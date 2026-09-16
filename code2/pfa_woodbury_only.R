source("basic_functions.R")

laplace_lambda_j_wbonly <- function(Y_j, X_j, M_j, mu, phi, B, sigma2,
                                    BtB, M_K, M_K_inv, log_det_MK,
                                    max_iter = 100, estep_gtol = 1e-3, bt_max = 30,
                                    c1 = 1e-4) {
  N_j <- nrow(Y_j)
  Q <- ncol(Y_j)
  K <- ncol(B)
  fixed <- matrix(rep(mu, each = N_j), N_j, Q) + X_j %*% phi
  Sinv <- function(v) v / sigma2 - B %*% (M_K_inv %*% (t(B) %*% v)) / sigma2^2
  lp <- function(a) {
    eta <- sweep(fixed, 2, a, "+")
    sum(Y_j * eta) - sum(M_j * row_logsumexp(eta)) - 0.5 * sum(a * Sinv(a))
  }

  wb_solve_exact <- function(g, pi_mat, d) {
    W <- sqrt(M_j) * pi_mat
    dt <- d + 1 / sigma2
    idt <- 1 / dt
    Widt <- sweep(W, 2, idt, "*")
    U <- cbind(t(W), B)
    Delta_11 <- diag(N_j) - Widt %*% t(W)
    Delta_12 <- -Widt %*% B
    Delta_22 <- sigma2^2 * M_K - crossprod(B, idt * B)
    Delta <- rbind(cbind(Delta_11, Delta_12), cbind(t(Delta_12), Delta_22))
    L <- tryCatch(chol(Delta + 1e-10 * diag(nrow(Delta))), error = function(e) NULL)
    if (is.null(L)) return(g * 0.01)
    idtg <- idt * g
    rhs <- as.numeric(crossprod(U, idtg))
    sol <- backsolve(L, forwardsolve(t(L), rhs))
    idtg + idt * as.numeric(U %*% sol)
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
    if (max(abs(g)) < estep_gtol) break
    dir <- wb_solve_exact(g, pi, d)
    gd <- sum(g * dir)
    step <- 1; accepted <- FALSE
    for (bt in 1:bt_max) {
      if (lp(lambda + step * dir) >= lp_cur + c1 * step * gd) { accepted <- TRUE; break }
      step <- step * 0.5
    }
    if (!accepted) {
      warning(sprintf("laplace_lambda_j_wbonly: line search failed to find an accepted step at inner iter %d (grad_norm=%.4g); returning current unconverged lambda", iter, max(abs(g))))
      break
    }
    lambda  <- lambda + step * dir
    lp_cur <- lp(lambda)
  }

  eta <- sweep(fixed, 2, lambda, "+")
  pi <- row_softmax(eta)
  Mpi <- sweep(pi, 1, M_j, "*")
  d <- as.numeric(colSums(Mpi))
  g <- as.numeric(colSums(Y_j - Mpi)) - Sinv(lambda)
  W <- sqrt(M_j) * pi
  dt <- d + 1 / sigma2
  idt <- 1 / dt
  Widt <- sweep(W, 2, idt, "*")
  U <- cbind(t(W), B)
  Delta_11 <- diag(N_j) - Widt %*% t(W)
  Delta_12 <- -Widt %*% B
  Delta_22 <- sigma2^2 * M_K - crossprod(B, idt * B)
  Delta <- rbind(cbind(Delta_11, Delta_12), cbind(t(Delta_12), Delta_22))
  L <- chol(Delta + 1e-10 * diag(nrow(Delta)))
  log_det_Delta <- 2 * sum(log(diag(L)))
  log_det_shat <- -sum(log(dt)) + log_det_MK + 2 * K * log(sigma2) - log_det_Delta
  DeltaInv <- chol2inv(L)
  UdiagIdt <- idt * U
  S_hat <- diag(idt, Q) + UdiagIdt %*% DeltaInv %*% t(UdiagIdt)

  list(lambda_hat = lambda, S_hat = S_hat, lp_mode = lp_cur,
       log_det_shat = log_det_shat, n_iter = iter, grad_norm = max(abs(g)))
}

estep_wbonly <- function(J, group, Y, X, M, mu, phi, B, sigma2,
                         max_iter = 100, estep_gtol = 1e-3) {
  Q <- length(mu)
  K <- ncol(B)
  BtB <- crossprod(B)
  M_K <- diag(K) + BtB / sigma2
  M_K_inv <- solve(M_K)
  log_det_MK <- 2 * sum(log(diag(chol(M_K))))
  ldSigma <- Q * log(sigma2) + log_det_MK

  lambda_hat <- matrix(0, Q, J); S_hat <- vector("list", J)
  lp_total <- 0; ld_S_total <- 0
  n_iter_vec <- numeric(J); grad_norm_vec <- numeric(J)

  for (j in 1:J) {
    idx <- which(group == j)
    if (!length(idx)) {
      S_hat[[j]] <- B %*% t(B) + sigma2 * diag(Q)
      ld_S_total <- ld_S_total + ldSigma
      n_iter_vec[j] <- NA; grad_norm_vec[j] <- NA
      next
    }

    res <- laplace_lambda_j_wbonly(Y[idx, , drop = FALSE], X[idx, , drop = FALSE],
                                   M[idx], mu, phi, B, sigma2,
                                   BtB = BtB, M_K = M_K, M_K_inv = M_K_inv,
                                   log_det_MK = log_det_MK,
                                   max_iter = max_iter, estep_gtol = estep_gtol)
    lambda_hat[, j] <- res$lambda_hat
    S_hat[[j]] <- res$S_hat
    lp_total <- lp_total + res$lp_mode
    ld_S_total <- ld_S_total + res$log_det_shat
    n_iter_vec[j] <- res$n_iter
    grad_norm_vec[j] <- res$grad_norm
  }

  list(lambda_hat = lambda_hat, S_hat = S_hat, lp_total = lp_total,
       ld_S_total = ld_S_total, log_det_Sigma = ldSigma,
       n_iter_vec = n_iter_vec, grad_norm_vec = grad_norm_vec)
}

mstep_phi_wbonly <- function(Y, X, group, lambda_hat, mu, phi, lambda_phi = 0) {
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

fit_pfa_wbonly <- function(Y, X, group, K,
                          M = rowSums(Y),
                          max_iter = 60, tol = 1e-4,
                          lambda_phi = 0, sigma2_init = 0.3,
                          verbose = FALSE,
                          estep_max_iter = 100, estep_gtol = 1e-3,
                          trace = FALSE, B_true = NULL) {

  N <- nrow(Y)
  Q <- ncol(Y)
  P <- ncol(X)
  J <- max(group)
  stopifnot(all(X[, 1] == 1))

  avg_prop <- colMeans(Y / pmax(rowSums(Y), 1))
  mu <- log(avg_prop + 1e-8)
  mu <- mu - mean(mu)
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
    u_frac_trace <- numeric(max_iter)
    best_iter <- NA_integer_
    best_log_ev <- -Inf
    best_state <- NULL
  }

  converged <- FALSE
  stop_reason <- "max_iter"
  em <- 0

  for (em in 1:max_iter) {

    es <- estep_wbonly(J, group, Y, X, M, mu, phi, B, sigma2,
                       max_iter = estep_max_iter, estep_gtol = estep_gtol)
    lambda_hat <- es$lambda_hat
    S_hat <- es$S_hat

    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total

    if (trace) monotone_trace[em] <- if (em == 1) TRUE else (log_ev[em] >= log_ev[em - 1])

    mp <- mstep_phi_wbonly(Y, X, group, lambda_hat, mu, phi, lambda_phi = lambda_phi)
    mu <- mp$mu
    phi <- mp$phi

    S_obs <- tcrossprod(lambda_hat) / J
    for (j in 1:J) S_obs <- S_obs + S_hat[[j]] / J
    rt <- ppca_closed(S_obs, K)
    B <- apply_PLT(rt$B)
    sigma2 <- rt$sigma2

    if (trace) {
      sigma2_trace[em] <- sigma2
      if (!is.null(B_true)) B_dist_trace[em] <- subspace_distance(B, B_true)
      u_frac_trace[em] <- u_frac(B)
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

    if (em > 1 && abs(log_ev[em] - log_ev[em - 1]) < tol) { converged <- TRUE; stop_reason <- "tol"; break }
  }

  es_final <- estep_wbonly(J, group, Y, X, M, mu, phi, B, sigma2,
                           max_iter = estep_max_iter, estep_gtol = estep_gtol)
  lambda_hat <- es_final$lambda_hat
  S_hat <- es_final$S_hat
  log_ev_final <- es_final$lp_total - 0.5 * J * es_final$log_det_Sigma + 0.5 * es_final$ld_S_total

  out <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
             lambda_hat = lambda_hat, S_hat = S_hat,
             log_evidence = log_ev[1:em], log_ev_final = log_ev_final,
             converged = converged, iterations = em, stop_reason = stop_reason,
             inner_iter_mean_final = mean(es_final$n_iter_vec, na.rm = TRUE),
             grad_norm_mean_final = mean(es_final$grad_norm_vec, na.rm = TRUE),
             inner_maxiter_hit_pct_final = mean(es_final$n_iter_vec >= estep_max_iter, na.rm = TRUE))

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
      u_frac_trace = u_frac_trace[1:em_final],
      best_iter = best_iter,
      best_log_ev = best_log_ev,
      best = best_state
    ))
  }

  out
}
