library(glmmTMB)
library(bench)
library(nnet) # multinomial
library(psych) # fa package
library(Matrix)

# glmmTMB reference fit (without idiosyncratic variance)
add_offset <- function(long, Y) {
  tot <- rowSums(Y)                  
  long$log_total <- log(pmax(tot[as.integer(long$obs)], 1))
  long
}

fit_pfa_glmmTMB <- function(Y, X, group, K,
                            covariate_terms = NULL, 
                            reml = FALSE) {
  long <- add_offset(build_long(Y, X, group), Y)

  rhs_fixed <- "category"
  if (!is.null(covariate_terms)) rhs_fixed <- paste(rhs_fixed, "+", covariate_terms)
  
  form <- as.formula(
    paste0("count ~ ", rhs_fixed,
           " + rr(category + 0 | group, d = ", K, ")")
  )
  
  glmmTMB(form,
          data = long,
          family = poisson(link = "log"),
          offset = long$log_total,
          REML = reml)
}

extract_loadings <- function(fit) {
  vc <- VarCorr(fit)$cond$group # Q x Q reduced-rank covariance factor loadings:
  # leading-K eigen-decomposition of the rr covariance
  e  <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
  # eigs_sym accelerate?
  K  <- attr(vc, "d"); if (is.null(K)) K <- sum(e$values > 1e-8)
  L  <- e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(e$values[1:K], 0)), K)
  list(B = L, cov = vc)
}

# dense multinomial follow the thesis
row_softmax <- function(eta) { mx <- apply(eta, 1, max); ee <- exp(eta - mx); ee / rowSums(ee) }
row_logsumexp <- function(eta) { mx <- apply(eta, 1, max); mx + log(rowSums(exp(eta - mx))) }

laplace_lambda_j_dense <- function(Y_j, X_j, M_j, mu, phi, Sigma_inv,
                                  max_iter = 100, gtol = 1e-3, bt_max = 30,
                                  c1 = 1e-4) {
  N_j <- nrow(Y_j); Q <- ncol(Y_j)
  fixed <- matrix(rep(mu, each = N_j), N_j, Q) + X_j %*% phi

  lp <- function(a) {
    eta <- sweep(fixed, 2, a, "+")
    sum(Y_j * eta) - sum(M_j * row_logsumexp(eta)) -
      0.5 * sum(a * (Sigma_inv %*% a))
  }

  lambda <- rep(0, Q)
  lp_cur <- lp(lambda)
  g <- rep(Inf, Q)

  for (iter in 1:max_iter) {
    eta <- sweep(fixed, 2, lambda, "+")
    pi <- row_softmax(eta)
    Mpi <- sweep(pi, 1, M_j, "*")
    g <- as.numeric(colSums(Y_j - Mpi)) - as.numeric(Sigma_inv %*% lambda)
    if (max(abs(g)) < gtol) break

    H <- diag(as.numeric(colSums(Mpi)), Q) - crossprod(pi, Mpi) + Sigma_inv
    dir <- tryCatch(solve(H, g), error = function(e) g * 0.01)  
    gd <- sum(g * dir)                                        

    step <- 1; accepted <- FALSE
    for (bt in 1:bt_max) {
      if (lp(lambda + step * dir) >= lp_cur + c1 * step * gd) { accepted <- TRUE; break }
      step <- step * 0.5
    }
    if (!accepted) break               
    lambda <- lambda + step * dir
    lp_cur <- lp(lambda)
  }

  eta <- sweep(fixed, 2, lambda, "+")
  pi <- row_softmax(eta)
  Mpi <- sweep(pi, 1, M_j, "*")
  H <- diag(as.numeric(colSums(Mpi)), Q) - crossprod(pi, Mpi) + Sigma_inv
  S_hat <- tryCatch({
    Hc <- Cholesky(H)
    S_hat <- solve(Hc, Diagonal(nrow(H)))
  }, error = function(e) {diag(Q) * 1e-3})
  ldS <- -as.numeric(determinant(H, logarithm = TRUE)$modulus)

  list(lambda_hat = lambda, S_hat = S_hat, lp_mode = lp_cur,
       log_det_shat = ldS, n_iter = iter, grad_norm = max(abs(g)))
}

estep_dense <- function(J, group, Y, X, M, mu, phi, B, sigma2,
                        max_iter = 100, gtol = 1e-3) {
  Q <- length(mu)
  Sigma <- B %*% t(B) + sigma2 * diag(Q)
  Sigma_inv <- tryCatch({
    cholSigma <- Cholesky(Sigma)
    solve(choSigma, Diagonal(Q))
  }, error = function(e) {diag(Q) / sigma2})
  ldSigma <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)

  lambda_hat <- matrix(0, Q, J)
  S_hat <- vector("list", J)
  lp_total <- 0
  ld_S_total <- 0

  for (j in 1:J) {
    idx <- which(group == j)
    if (!length(idx)) {
      S_hat[[j]] <- Sigma
      ld_S_total <- ld_S_total + ldSigma
      next
    }
    res <- laplace_lambda_j_dense(Y[idx, , drop = FALSE], X[idx, , drop = FALSE],
                                 M[idx], mu, phi, Sigma_inv,
                                 max_iter = max_iter, gtol = gtol)
    lambda_hat[, j] <- res$lambda_hat
    S_hat[[j]] <- res$S_hat
    lp_total <- lp_total + res$lp_mode
    ld_S_total <- ld_S_total + res$log_det_shat
  }
  list(lambda_hat = lambda_hat, S_hat = S_hat,
       lp_total = lp_total, ld_S_total = ld_S_total,
       log_det_Sigma = ldSigma)
}

mstep_phi_dense <- function(Y, X, group, lambda_hat, mu, phi, lambda_phi = 0) {
  N <- nrow(Y); Q <- ncol(Y); P <- ncol(X)
  A_obs <- t(lambda_hat[, group, drop = FALSE])        

  eta <- sweep(X %*% phi, 2, mu, "+") + A_obs
  M_i <- rowSums(Y)
  delta <- log(pmax(M_i, 1)) - row_logsumexp(eta)

  mu_new <- mu
  phi_new <- phi

  for (q in 1:Q) {
    off <- delta + A_obs[, q]
    co <- tryCatch({
      if (lambda_phi <= 0) {
        fit <- suppressWarnings(
          glm.fit(x = X, y = Y[, q], family = poisson(), offset = off))
        fit$coefficients
      } else {
        pois_ridge_irls(X, Y[, q], off, lambda_phi)
      }
    }, error = function(e) c(mu[q], phi[-1, q]))
    if (any(!is.finite(co))) co <- c(mu[q], phi[-1, q])
    mu_new[q] <- co[1]
    if (P > 1) phi_new[2:P, q] <- co[2:P]
  }
  phi_new[1, ] <- 0
  list(mu = mu_new, phi = phi_new)
}

pois_ridge_irls <- function(X, y, off, lambda, max_iter = 50, tol = 1e-8) {
  P <- ncol(X)
  pen <- diag(c(0, rep(lambda, P - 1)), P)
  beta <- rep(0, P)
  for (it in 1:max_iter) {
    eta <- as.numeric(X %*% beta) + off
    mu_ <- exp(pmin(eta, 30))
    W <- mu_
    z <- eta - off + (y - mu_) / pmax(mu_, 1e-10)
    XtW <- t(X * W)
    b_new <- tryCatch(solve(XtW %*% X + pen, XtW %*% z),
                      error = function(e) beta)
    if (max(abs(b_new - beta)) < tol) { beta <- b_new; break }
    beta <- b_new
  }
  as.numeric(beta)
}

rubin_thayer_spherical <- function(Sigma_obs, K, B_init = NULL,
                                   sigma2_init = 0.3,
                                   max_iter = 500, tol = 1e-10) {
  Q <- nrow(Sigma_obs)
  if (is.null(B_init)) {
    sv <- svd(Sigma_obs, nu = K, nv = K)
    lam <- pmax(sv$d[1:K] - sigma2_init, 0.05)
    B <- sv$u %*% diag(sqrt(lam), K)
  } else B <- B_init
  sigma2 <- max(sigma2_init, 1e-6)

  for (iter in 1:max_iter) {
    B_old <- B; s_old <- sigma2
    Sigma <- B %*% t(B) + sigma2 * diag(Q)
    Si <- tryCatch(solve(Sigma), error = function(e) diag(Q) / sigma2)
    beta <- t(B) %*% Si                                   # K x Q
    Theta <- diag(K) - beta %*% B + beta %*% Sigma_obs %*% t(beta)
    B_new <- Sigma_obs %*% t(beta) %*% solve(Theta)
    sigma2_new <- max(sum(diag(Sigma_obs - B_new %*% beta %*% Sigma_obs)) / Q,
                      1e-6)
    B <- B_new; sigma2 <- sigma2_new
    if (max(abs(B - B_old)) < tol && abs(sigma2 - s_old) < tol) break
  }
  list(B = B, sigma2 = sigma2)
}

fit_pfa_dense <- function(Y, X, group, K,
                          M = rowSums(Y),
                          max_iter = 60, tol = 1e-4,
                          lambda_phi = 0, sigma2_init = 0.3,
                          verbose = FALSE,
                          estep_max_iter = 100, estep_gtol = 1e-3) {

  N <- nrow(Y); Q <- ncol(Y); P <- ncol(X); J <- max(group)
  stopifnot(all(X[, 1] == 1))

  avg_prop <- colMeans(Y / pmax(rowSums(Y), 1))
  mu <- log(avg_prop + 1e-8); mu <- mu - mean(mu)
  phi <- matrix(0, P, Q)

  gm <- matrix(0, J, Q)                         
  for (j in 1:J) {
    idx <- which(group == j)
    if (length(idx) > 0) {
      gs <- colSums(Y[idx, , drop = FALSE])
      gm[j, ] <- log((gs + 1e-5) / (sum(gs) + Q * 1e-5))
    }
  }
  gm_c <- sweep(gm, 2, colMeans(gm))                 
  sv0 <- svd(t(gm_c), nu = K, nv = K)          
  B <- sv0$u %*% diag(pmax(sv0$d[1:K] * 0.5, 0.1), K)
  B <- apply_PLT(B)
  sigma2 <- sigma2_init

  log_ev <- numeric(max_iter)
  converged <- FALSE
  em <- 0

  for (em in 1:max_iter) {

    es <- estep_dense(J, group, Y, X, M, mu, phi, B, sigma2,
                             max_iter = estep_max_iter, gtol = estep_gtol)
    lambda_hat <- es$lambda_hat
    S_hat <- es$S_hat

    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma +
      0.5 * es$ld_S_total

    mp <- mstep_phi_dense(Y, X, group, lambda_hat, mu, phi,
                           lambda_phi = lambda_phi)
    mu <- mp$mu; phi <- mp$phi

    S_obs <- tcrossprod(lambda_hat) / J
    for (j in 1:J) S_obs <- S_obs + S_hat[[j]] / J
    rt <- rubin_thayer_spherical(S_obs, K,
                                     B_init = B, sigma2_init = sigma2)
    B <- apply_PLT(rt$B)
    sigma2 <- rt$sigma2

    if (verbose) cat(sprintf("iter %3d  log_ev = %.4f  sigma2 = %.4f\n",
                             em, log_ev[em], sigma2))

    if (em > 1 &&
        abs(log_ev[em] - log_ev[em - 1]) <
        tol * (abs(log_ev[em - 1]) + 1)) { converged <- TRUE; break }
  }

  list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
       lambda_hat = lambda_hat, S_hat = S_hat,
       log_evidence = log_ev[1:em], converged = converged, iterations = em)
}
