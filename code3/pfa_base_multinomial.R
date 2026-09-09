laplace_lambda_j_dense <- function(Y_j, X_j, M_j, mu, phi, Sigma_inv,
                                   max_iter = 100, estep_gtol = 1e-3, bt_max = 30,
                                   c1 = 1e-4) {
  N_j <- nrow(Y_j)
  Q <- ncol(Y_j)
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
    if (max(abs(g)) < estep_gtol) break
    
    H <- diag(as.numeric(colSums(Mpi)), Q) - crossprod(pi, Mpi) + Sigma_inv
    dir <- tryCatch(solve(H, g), error = function(e) g * 0.01)
    gd <- sum(g * dir)
    
    step <- 1
    accepted <- FALSE
    for (bt in 1:bt_max) {
      if (lp(lambda + step * dir) >= lp_cur + c1 * step * gd) { accepted <- TRUE; break }
      step <- step * 0.5
    }
    if (!accepted) {
      warning(sprintf("laplace_lambda_j_dense: line search failed to find an accepted step at inner iter %d (grad_norm=%.4g); returning current unconverged lambda", iter, max(abs(g))))
      break
    }
    lambda <- lambda + step * dir
    lp_cur <- lp(lambda)
  }
  
  eta <- sweep(fixed, 2, lambda, "+")
  pi <- row_softmax(eta)
  Mpi <- sweep(pi, 1, M_j, "*")
  H <- diag(as.numeric(colSums(Mpi)), Q) - crossprod(pi, Mpi) + Sigma_inv
  S_hat <- tryCatch({
    Hc <- Cholesky(H)
    solve(Hc, Diagonal(nrow(H)))
  }, error = function(e) {
    warning(sprintf("Cholesky (H) error: %s", conditionMessage(e)))
    diag(Q) * 1e-3
  })
  ldS <- -as.numeric(determinant(H, logarithm = TRUE)$modulus)
  
  list(lambda_hat = lambda, S_hat = S_hat, lp_mode = lp_cur,
       log_det_shat = ldS, n_iter = iter, grad_norm = max(abs(g)))
}

estep_dense <- function(J, group, Y, X, M, mu, phi, B, sigma2,
                        max_iter = 100, estep_gtol = 1e-3) {
  Q <- length(mu)
  
  Sigma <- B %*% t(B) + sigma2 * diag(Q)
  Sigma_inv <- tryCatch({
    cholSigma <- Cholesky(Sigma)
    solve(cholSigma, Diagonal(Q))
  }, error = function(e) {
    warning(sprintf("Cholesky (Sigma) error for sigma2=%.6g: %s", sigma2, conditionMessage(e)))
    diag(Q) / sigma2})
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
                                  max_iter = max_iter, estep_gtol = estep_gtol)
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

fit_pfa_dense <- function(Y, X, group, K,
                          M = rowSums(Y),
                          max_iter = 60, tol = 1e-4,
                          lambda_phi = 0, sigma2_init = 0.3,
                          verbose = FALSE,
                          estep_max_iter = 100, estep_gtol = 1e-3,
                          trace = FALSE, B_true = NULL) {
  
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
  if (trace) {
    sigma2_trace <- numeric(max_iter)
    monotone_trace <- rep(NA, max_iter)
    B_dist_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
    best_iter <- NA_integer_
    best_log_ev <- -Inf
    best_state <- NULL
  }
  
  converged <- FALSE
  stop_reason <- "max_iter"
  em <- 0
  
  for (em in 1:max_iter) {
    
    es <- estep_dense(J, group, Y, X, M, mu, phi, B, sigma2,
                      max_iter = estep_max_iter, estep_gtol = estep_gtol)
    lambda_hat <- es$lambda_hat
    S_hat <- es$S_hat
    
    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma +
      0.5 * es$ld_S_total
    
    if (trace) monotone_trace[em] <- if (em == 1) TRUE else (log_ev[em] >= log_ev[em - 1])
    
    mp <- mstep_phi_dense(Y, X, group, lambda_hat, mu, phi,
                          lambda_phi = lambda_phi)
    mu <- mp$mu
    phi <- mp$phi
    
    S_obs <- tcrossprod(lambda_hat) / J
    for (j in 1:J) S_obs <- S_obs + S_hat[[j]] / J
    rt <- ppca(S_obs, K)
    B <- apply_PLT(rt$B)
    sigma2 <- rt$sigma2
    
    if (trace) {
      sigma2_trace[em] <- sigma2
      if (!is.null(B_true)) B_dist_trace[em] <- subspace_dist(B, B_true)
      if (log_ev[em] > best_log_ev) {
        best_log_ev <- log_ev[em]
        best_iter <- em
        best_state <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
                           lambda_hat = lambda_hat, S_hat = S_hat)
      }
    }
    
    if (verbose) {
      bd_str <- if (trace && !is.null(B_true)) sprintf("  B_dist = %.4f", B_dist_trace[em]) else ""
      cat(sprintf("iter %3d  log_ev = %.4f  sigma2 = %.4f%s\n",
                  em, log_ev[em], sigma2, bd_str))
    }
    
    if (em > 1 &&
        abs(log_ev[em] - log_ev[em - 1]) < tol) { converged <- TRUE; stop_reason <- "tol"; break }
  }
  
  es_final <- estep_dense(J, group, Y, X, M, mu, phi, B, sigma2,
                          max_iter = estep_max_iter, estep_gtol = estep_gtol)
  lambda_hat <- es_final$lambda_hat
  S_hat <- es_final$S_hat
  log_ev_final <- es_final$lp_total - 0.5 * J * es_final$log_det_Sigma +
    0.5 * es_final$ld_S_total
  
  out <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
              lambda_hat = lambda_hat, S_hat = S_hat,
              log_evidence = log_ev[1:em], log_ev_final = log_ev_final,
              converged = converged, iterations = em, stop_reason = stop_reason)
  
  if (trace) {
    em_final <- em
    is_monotone <- all(monotone_trace[1:em_final])
    first_decrease_iter <- if (!is_monotone) which(!monotone_trace[1:em_final])[1] else NA_integer_
    out <- c(out, list(
      sigma2_trace = sigma2_trace[1:em_final],
      monotone_trace = monotone_trace[1:em_final],
      is_monotone = is_monotone,
      first_decrease_iter = first_decrease_iter,
      B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em_final] else NULL,
      best_iter = best_iter,
      best_log_ev = best_log_ev,
      best = best_state
    ))
  }
  
  out
}