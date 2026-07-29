# =============================================================================
# pfa_woodbury_armijo.R
# =============================================================================


# =============================================================================
# 1. SHARED HELPERS  (identical to the dense file)
# =============================================================================

apply_PLT <- function(B) {
  K <- ncol(B)
  qr_obj <- qr(t(B[1:K, , drop = FALSE]))
  B_rot <- B %*% qr.Q(qr_obj)
  for (k in 1:K) if (B_rot[k, k] < 0) B_rot[, k] <- -B_rot[, k]
  for (k in 1:K) if (k > 1) B_rot[1:(k - 1), k] <- 0
  B_rot
}

row_softmax <- function(eta) { mx <- apply(eta, 1, max); ee <- exp(eta - mx); ee / rowSums(ee) }
row_logsumexp <- function(eta) { mx <- apply(eta, 1, max); mx + log(rowSums(exp(eta - mx))) }


# =============================================================================
# 2. E-STEP: Woodbury Laplace update for one group
# =============================================================================

laplace_lambda_j_wb2 <- function(Y_j, X_j, M_j, mu, phi, B, sigma2,
                                BtB, M_K, M_K_inv, log_det_MK,
                                max_iter = 100, gtol = 1e-3, bt_max = 30,
                                c1 = 1e-4, exact_Shat = FALSE) {
  N_j <- nrow(Y_j); Q <- ncol(Y_j); K <- ncol(B)
  fixed <- matrix(rep(mu, each = N_j), N_j, Q) + X_j %*% phi

  Sinv <- function(v) v / sigma2 - B %*% (M_K_inv %*% (t(B) %*% v)) / sigma2^2

  lp <- function(a) {
    eta <- sweep(fixed, 2, a, "+")
    sum(Y_j * eta) - sum(M_j * row_logsumexp(eta)) - 0.5 * sum(a * Sinv(a))
  }

  wb_solve <- function(g, d) {
    dt <- d + 1 / sigma2             
    idt <- 1 / dt
    Gam <- t(B) %*% (idt * B) - sigma2^2 * M_K
    L <- tryCatch(chol(-Gam + 1e-10 * diag(K)), error = function(e) NULL)
    if (is.null(L)) return(g * 0.01)  
    idt * g - idt * (B %*% backsolve(L, forwardsolve(t(L), -crossprod(B, idt * g))))
  }

  lambda <- rep(0, Q)
  lp_cur <- lp(lambda)
  g <- rep(Inf, Q)
  d <- rep(0, Q)

  for (iter in 1:max_iter) {
    eta <- sweep(fixed, 2, lambda, "+")
    pi <- row_softmax(eta)
    Mpi <- sweep(pi, 1, M_j, "*")

    d <- as.numeric(colSums(Mpi))                       
    g <- as.numeric(colSums(Y_j - Mpi)) - Sinv(lambda)   

    if (max(abs(g)) < gtol) break

    dir <- wb_solve(g, d)                           
    gd <- sum(g * dir)                           

    step <- 1; accepted <- FALSE
    for (bt in 1:bt_max) {
      if (lp(lambda + step * dir) >= lp_cur + c1 * step * gd) { accepted <- TRUE; break }
      step <- step * 0.5
    }
    if (!accepted) break             
    lambda  <- lambda + step * dir
    lp_cur <- lp(lambda)
  }

  if (exact_Shat) {
    Sinv_mat <- diag(1 / sigma2, Q) - B %*% M_K_inv %*% t(B) / sigma2^2
    eta <- sweep(fixed, 2, lambda, "+")
    pi  <- row_softmax(eta)
    Mpi <- sweep(pi, 1, M_j, "*")
    negH <- diag(as.numeric(colSums(Mpi)), Q) - crossprod(pi, Mpi) + Sinv_mat
    cL  <- chol(negH)                                  
    S_hat        <- chol2inv(cL)
    log_det_shat <- -2 * sum(log(diag(cL)))
  } else {
    dt  <- d + 1 / sigma2
    idt <- 1 / dt
    Gam <- t(B) %*% (idt * B) - sigma2^2 * M_K
    L   <- chol(-Gam + 1e-10 * diag(K))
    log_det_shat <- -(sum(log(dt)) - 2 * K * log(sigma2) - log_det_MK +
                        2 * sum(log(diag(L))))
    GamInv <- chol2inv(L)
    IdtB   <- idt * B
    S_hat  <- diag(idt, Q) + IdtB %*% GamInv %*% t(IdtB)
  }
  
  list(lambda_hat = lambda, S_hat = S_hat, lp_mode = lp_cur,
       log_det_shat = log_det_shat, n_iter = iter, grad_norm = max(abs(g)))
}

estep_wb <- function(J, group, Y, X, M, mu, phi, B, sigma2,
                     max_iter = 100, gtol = 1e-3, exact_Shat = FALSE) {
  Q <- length(mu); K <- ncol(B)

  BtB <- crossprod(B)
  M_K <- diag(K) + BtB / sigma2
  M_K_inv <- solve(M_K)
  log_det_MK <- 2 * sum(log(diag(chol(M_K))))
  ldSigma <- Q * log(sigma2) + log_det_MK         

  lambda_hat <- matrix(0, Q, J); S_hat <- vector("list", J)
  lp_total  <- 0; ld_S_total <- 0

  for (j in 1:J) {
    idx <- which(group == j)
    if (!length(idx)) {
      S_hat[[j]] <- B %*% t(B) + sigma2 * diag(Q)       
      ld_S_total <- ld_S_total + ldSigma
      next
    }
    res <- laplace_lambda_j_wb2(Y[idx, , drop = FALSE], X[idx, , drop = FALSE],
                               M[idx], mu, phi, B, sigma2,
                               BtB = BtB, M_K = M_K, M_K_inv = M_K_inv,
                               log_det_MK = log_det_MK,
                               max_iter = max_iter, gtol = gtol, exact_Shat = exact_Shat)
    lambda_hat[, j] <- res$lambda_hat
    S_hat[[j]] <- res$S_hat
    lp_total <- lp_total + res$lp_mode
    ld_S_total <- ld_S_total + res$log_det_shat
  }
  list(lambda_hat = lambda_hat, S_hat = S_hat, lp_total = lp_total,
       ld_S_total = ld_S_total, log_det_Sigma = ldSigma)
}


# =============================================================================
# 3. M-STEP (a): fixed effects (mu, phi) 
# =============================================================================

mstep_phi_wb <- function(Y, X, group, lambda_hat, mu, phi, lambda_phi = 0) {
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
      fit <- suppressWarnings(
        glm.fit(x = X, y = Y[, q], family = poisson(), offset = off))
      fit$coefficients
    }, error = function(e) c(mu[q], phi[-1, q]))
    if (any(!is.finite(co))) co <- c(mu[q], phi[-1, q])
    mu_new[q] <- co[1]
    if (P > 1) phi_new[2:P, q] <- co[2:P]
  }
  phi_new[1, ] <- 0
  list(mu = mu_new, phi = phi_new)
}


# =============================================================================
# 4. M-STEP (b): WOODBURY RUBIN-THAYER
# =============================================================================

rubin_thayer_wb <- function(Sigma_obs, K, B_init = NULL, sigma2_init = 0.3,
                            max_iter = 500, tol = 1e-10) {
  Q <- nrow(Sigma_obs)
  trS <- sum(diag(Sigma_obs))

  if (is.null(B_init)) {
    sv <- svd(Sigma_obs, nu = K, nv = K)
    lam <- pmax(sv$d[1:K] - sigma2_init, 0.05)
    B <- sv$u %*% diag(sqrt(lam), K)
  } else B <- B_init
  sigma2 <- max(sigma2_init, 1e-6)

  for (iter in 1:max_iter) {
    B_old <- B; s_old <- sigma2
    BtB_ <- crossprod(B)
    M_K_ <- diag(K) + BtB_ / sigma2
    Mki <- solve(M_K_)                     

    SoB <- Sigma_obs %*% B             
    BtSB <- crossprod(B, SoB)     

    # Theta = M_K^{-1} + (1/sigma2^2) * M_K^{-1} * BtSB * M_K^{-1}
    Theta  <- Mki + Mki %*% BtSB %*% Mki / sigma2^2   # K x K
    B_new  <- SoB %*% (Mki %*% solve(Theta)) / sigma2  # Q x K

    # sigma2 = (tr(Sigma_obs) - (1/sigma2)*tr(M_K^{-1} B' Sigma_obs B_new)) / Q
    BtSBn <- crossprod(SoB, B_new)           # K x K: B' Sigma_obs B_new
    trace_term <- sum(Mki * t(BtSBn)) / sigma2
    sigma2_new <- max((trS - trace_term) / Q, 1e-6)

    B <- B_new; sigma2 <- sigma2_new
    if (max(abs(B - B_old)) < tol && abs(sigma2 - s_old) < tol) break
  }
  list(B = B, sigma2 = sigma2)
}


# =============================================================================
# 5. MAIN FITTER  (fully Woodbury-accelerated)
# =============================================================================

fit_pfa_woodbury <- function(Y, X, group, K,
                             M = rowSums(Y),
                             max_iter = 60, tol = 1e-4,
                             lambda_phi = 0, sigma2_init = 0.3,
                             verbose = FALSE,
                             estep_max_iter = 100, estep_gtol = 1e-3, exact_Shat = FALSE) {

  N <- nrow(Y); Q <- ncol(Y); P <- ncol(X); J <- max(group)
  stopifnot(all(X[, 1] == 1))

  # ------------------------------------------------------------------
  # Initialisation
  # ------------------------------------------------------------------
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

    # ---- E-step ----
    es <- estep_wb(J, group, Y, X, M, mu, phi, B, sigma2,
                          max_iter = estep_max_iter, gtol = estep_gtol, exact_Shat = exact_Shat)
    lambda_hat <- es$lambda_hat
    S_hat <- es$S_hat

    # ---- Laplace log-evidence monitor ----
    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma +
      0.5 * es$ld_S_total

    # ---- M-step (a): fixed effects ----
    mp <- mstep_phi_wb(Y, X, group, lambda_hat, mu, phi,
                        lambda_phi = lambda_phi)
    mu <- mp$mu; phi <- mp$phi

    # ---- M-step (b): factor loadings ----
    S_obs <- tcrossprod(lambda_hat) / J
    for (j in 1:J) S_obs <- S_obs + S_hat[[j]] / J
    rt <- rubin_thayer_wb(S_obs, K, B_init = B, sigma2_init = sigma2)
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
