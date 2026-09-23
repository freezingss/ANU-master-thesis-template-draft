source("setup.R")
source("basic_functions.R")

laplace_lambda_j_wbonly_warm_start <- function(Y_j, X_j, M_j, mu, phi, B, sigma2,
                                    BtB, M_K, M_K_inv, log_det_MK,
                                    lambda_init = NULL,
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

  lambda <- if (is.null(lambda_init)) rep(0, Q) else lambda_init
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
      warning(sprintf("laplace_lambda_j_wbonly_warm_start: line search failed to find an accepted step at inner iter %d (grad_norm=%.4g); returning current unconverged lambda", iter, max(abs(g))))
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

estep_wbonly_warm_start <- function(J, group, Y, X, M, mu, phi, B, sigma2,
                                    lambda_init = NULL,
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

    res <- laplace_lambda_j_wbonly_warm_start(Y[idx, , drop = FALSE], X[idx, , drop = FALSE],
                                   M[idx], mu, phi, B, sigma2,
                                   BtB = BtB, M_K = M_K, M_K_inv = M_K_inv,
                                   log_det_MK = log_det_MK,
                                   lambda_init = if (is.null(lambda_init)) NULL else lambda_init[, j],
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

mstep_phi_wbonly_warm_start <- function(Y, X, group, lambda_hat, mu, phi, lambda_phi = 0) {
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

  mu_new <- mu_new - mean(mu_new)
  if (P > 1) phi_new[2:P, ] <- sweep(phi_new[2:P, , drop = FALSE], 1, rowMeans(phi_new[2:P, , drop = FALSE]))

  list(mu = mu_new, phi = phi_new)
}

extrapolate_remaining <- function(shift_hist, em, burn_in, rate_window) {
  if (em <= burn_in) return(list(r_hat = NA_real_, remaining = NA_real_))
  x <- pmax(shift_hist[1:em], 1e-12)
  n_x <- length(x)
  idx <- max(1, n_x - rate_window + 1):n_x
  xw <- x[idx]
  rate_fit <- tryCatch(lm(log(xw) ~ idx), error = function(e) NULL)
  if (is.null(rate_fit)) return(list(r_hat = NA_real_, remaining = NA_real_))
  r_hat <- exp(coef(rate_fit)[2])
  if (!is.finite(r_hat) || r_hat <= 0 || r_hat >= 1) return(list(r_hat = NA_real_, remaining = NA_real_))
  remaining <- xw[length(xw)] * r_hat / (1 - r_hat)
  list(r_hat = r_hat, remaining = remaining)
}

default_tols <- function(mu_init, sigma2_init, K, rel_tol = 1e-3) {
  mu_scale <- max(mu_init) - min(mu_init)
  if (!is.finite(mu_scale) || mu_scale <= 0) mu_scale <- 1
  list(tol_mu = rel_tol * mu_scale,
      tol_phi = rel_tol * mu_scale,
      tol_sigma2 = rel_tol * sigma2_init,
      tol_B = rel_tol * sqrt(K))
}

fit_pfa_wbonly_warm_start <- function(Y, X, group, K,
                                      M = rowSums(Y),
                                      max_iter = 60, tol = 1e-4,
                                      tol_mu = NULL, tol_phi = NULL,
                                      tol_sigma2 = NULL, tol_B = NULL,
                                      rel_tol = 1e-3,
                                      lambda_phi = 0, sigma2_init = 0.3,
                                      verbose = FALSE,
                                      estep_max_iter = 100, estep_gtol = 1e-3,
                                      trace = FALSE, B_true = NULL,
                                      patience = 15, burn_in = 15, rate_window = 10,
                                      stall_patience = 5, also_require_mu_phi = FALSE) {

  N <- nrow(Y)
  Q <- ncol(Y)
  P <- ncol(X)
  J <- max(group)
  stopifnot(all(X[, 1] == 1))

  avg_prop <- colMeans(Y / pmax(rowSums(Y), 1))
  mu <- log(avg_prop + 1e-8)
  mu <- mu - mean(mu)
  phi <- matrix(0, P, Q)

  dt0 <- default_tols(mu, sigma2_init, K, rel_tol = rel_tol)
  if (is.null(tol_mu)) tol_mu <- dt0$tol_mu
  if (is.null(tol_phi)) tol_phi <- dt0$tol_phi
  if (is.null(tol_sigma2)) tol_sigma2 <- dt0$tol_sigma2
  if (is.null(tol_B)) tol_B <- dt0$tol_B

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

  lambda_hat <- matrix(0, Q, J)

  log_ev <- numeric(max_iter)
  mu_shift_hist <- numeric(max_iter)
  phi_shift_hist <- numeric(max_iter)
  B_shift_hist <- numeric(max_iter)
  sigma2_shift_hist <- numeric(max_iter)

  if (trace) {
    mu_bar_trace <- numeric(max_iter)
    phi_bar_trace <- matrix(0, P - 1, max_iter)
    sigma2_trace <- numeric(max_iter)
    monotone_trace <- rep(NA, max_iter)
    B_dist_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
    u_frac_trace <- numeric(max_iter)
    r_hat_trace <- rep(NA_real_, max_iter)
    remaining_gap_trace <- rep(NA_real_, max_iter)
    inner_iter_mean_trace <- rep(NA_real_, max_iter)
    mu_r_hat_trace <- rep(NA_real_, max_iter)
    phi_r_hat_trace <- rep(NA_real_, max_iter)
    B_r_hat_trace <- rep(NA_real_, max_iter)
    sigma2_r_hat_trace <- rep(NA_real_, max_iter)
    mu_remaining_trace <- rep(NA_real_, max_iter)
    phi_remaining_trace <- rep(NA_real_, max_iter)
    B_remaining_trace <- rep(NA_real_, max_iter)
    sigma2_remaining_trace <- rep(NA_real_, max_iter)
    no_stall_trace <- rep(NA_integer_, max_iter)
    ev_ok_trace <- rep(FALSE, max_iter)
  }

  best_log_ev <- -Inf
  best_iter <- NA_integer_
  best_state <- NULL
  no_improve <- 0
  no_stall <- 0

  converged <- FALSE
  stop_reason <- "max_iter"
  em <- 0

  for (em in 1:max_iter) {

    es <- estep_wbonly_warm_start(J, group, Y, X, M, mu, phi, B, sigma2,
                                  lambda_init = lambda_hat,
                                  max_iter = estep_max_iter, estep_gtol = estep_gtol)
    lambda_hat <- es$lambda_hat
    S_hat <- es$S_hat
    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total

    if (trace) monotone_trace[em] <- if (em == 1) TRUE else (log_ev[em] >= log_ev[em - 1])

    if (log_ev[em] > best_log_ev) {
      best_log_ev <- log_ev[em]
      best_iter <- em
      best_state <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
                         lambda_hat = lambda_hat, S_hat = S_hat)
      no_improve <- 0
    } else {
      no_improve <- no_improve + 1
    }

    mu_prev <- mu
    phi_prev <- phi

    mp <- mstep_phi_wbonly_warm_start(Y, X, group, lambda_hat, mu, phi, lambda_phi = lambda_phi)
    mu <- mp$mu
    phi <- mp$phi

    mu_shift <- max(abs(mu - mu_prev))
    phi_shift <- if (P > 1) max(abs(phi[2:P, ] - phi_prev[2:P, , drop = FALSE])) else 0

    B_prev <- B
    sigma2_prev <- sigma2

    S_obs <- tcrossprod(lambda_hat) / J
    for (j in 1:J) S_obs <- S_obs + S_hat[[j]] / J
    rt <- ppca_closed(S_obs, K)
    B <- apply_PLT(rt$B)
    sigma2 <- rt$sigma2

    B_shift <- subspace_distance(B, B_prev)
    sigma2_shift <- abs(sigma2 - sigma2_prev)

    mu_shift_hist[em] <- mu_shift
    phi_shift_hist[em] <- phi_shift
    B_shift_hist[em] <- B_shift
    sigma2_shift_hist[em] <- sigma2_shift

    mu_ext <- extrapolate_remaining(mu_shift_hist, em, burn_in, rate_window)
    phi_ext <- extrapolate_remaining(phi_shift_hist, em, burn_in, rate_window)
    B_ext <- extrapolate_remaining(B_shift_hist, em, burn_in, rate_window)
    sigma2_ext <- extrapolate_remaining(sigma2_shift_hist, em, burn_in, rate_window)

    mu_ok <- if (!is.na(mu_ext$remaining)) mu_ext$remaining < tol_mu else mu_shift < tol_mu
    phi_ok <- if (!is.na(phi_ext$remaining)) phi_ext$remaining < tol_phi else phi_shift < tol_phi
    B_ok <- if (!is.na(B_ext$remaining)) B_ext$remaining < tol_B else B_shift < tol_B
    sigma2_ok <- if (!is.na(sigma2_ext$remaining)) sigma2_ext$remaining < tol_sigma2 else sigma2_shift < tol_sigma2

    if (B_ok && sigma2_ok && (!also_require_mu_phi || (mu_ok && phi_ok))) {
      no_stall <- no_stall + 1
    } else {
      no_stall <- 0
    }

    if (trace) {
      mu_bar_trace[em] <- mean(mu)
      if (P > 1) phi_bar_trace[, em] <- rowMeans(phi[2:P, , drop = FALSE])
      sigma2_trace[em] <- sigma2
      if (!is.null(B_true)) B_dist_trace[em] <- subspace_distance(B, B_true)
      u_frac_trace[em] <- u_frac(B)
      inner_iter_mean_trace[em] <- mean(es$n_iter_vec, na.rm = TRUE)
      mu_r_hat_trace[em] <- mu_ext$r_hat
      phi_r_hat_trace[em] <- phi_ext$r_hat
      B_r_hat_trace[em] <- B_ext$r_hat
      sigma2_r_hat_trace[em] <- sigma2_ext$r_hat
      mu_remaining_trace[em] <- mu_ext$remaining
      phi_remaining_trace[em] <- phi_ext$remaining
      B_remaining_trace[em] <- B_ext$remaining
      sigma2_remaining_trace[em] <- sigma2_ext$remaining
      no_stall_trace[em] <- no_stall
    }

    ev_ok <- FALSE
    r_hat <- NA_real_
    remaining_gap <- NA_real_
    if (em > burn_in) {
      d <- diff(log_ev[1:em])
      n_d <- length(d)
      idx <- max(1, n_d - rate_window + 1):n_d
      dw <- d[idx]
      if (all(dw > 0)) {
        rate_fit <- tryCatch(lm(log(dw) ~ idx), error = function(e) NULL)
        if (!is.null(rate_fit)) {
          r_hat <- exp(coef(rate_fit)[2])
          if (is.finite(r_hat) && r_hat > 0 && r_hat < 1) {
            remaining_gap <- dw[length(dw)] * r_hat / (1 - r_hat)
            ev_ok <- remaining_gap < tol
          }
        }
      }
    }

    if (trace) {
      r_hat_trace[em] <- r_hat
      remaining_gap_trace[em] <- remaining_gap
      ev_ok_trace[em] <- ev_ok
    }

    if (verbose) {
      cat(sprintf("iter %3d  B_shift = %.6f  sigma2_shift = %.6f  mu_shift = %.6f  phi_shift = %.6f  no_stall = %d\n",
                  em, B_shift, sigma2_shift, mu_shift, phi_shift, no_stall))
    }

    if (no_stall >= stall_patience) {
      converged <- TRUE; stop_reason <- "param_stable"; break
    }
    if (em > 1 && abs(log_ev[em] - log_ev[em - 1]) < tol) {
      converged <- TRUE; stop_reason <- "tol"; break
    }
    if (no_improve >= patience) {
      stop_reason <- "patience"; break
    }
  }

  out <- list(mu = best_state$mu, phi = best_state$phi,
             B = best_state$B, sigma2 = best_state$sigma2,
             lambda_hat = best_state$lambda_hat, S_hat = best_state$S_hat,
             log_evidence = log_ev[1:em], best_iter = best_iter, best_log_ev = best_log_ev,
             converged = converged, iterations = em, stop_reason = stop_reason,
             tol_mu = tol_mu, tol_phi = tol_phi, tol_sigma2 = tol_sigma2, tol_B = tol_B)

  if (trace) {
    em_final <- em
    is_monotone <- all(monotone_trace[1:em_final])
    first_decrease_iter <- if (!is_monotone) which(!monotone_trace[1:em_final])[1] else NA_integer_
    out <- c(out, list(
      mu_bar_trace = mu_bar_trace[1:em_final],
      phi_bar_trace = phi_bar_trace[, 1:em_final, drop = FALSE],
      sigma2_trace = sigma2_trace[1:em_final],
      monotone_trace = monotone_trace[1:em_final],
      is_monotone = is_monotone,
      first_decrease_iter = first_decrease_iter,
      B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em_final] else NULL,
      u_frac_trace = u_frac_trace[1:em_final],
      r_hat_trace = r_hat_trace[1:em_final],
      remaining_gap_trace = remaining_gap_trace[1:em_final],
      inner_iter_mean_trace = inner_iter_mean_trace[1:em_final],
      mu_shift_trace = mu_shift_hist[1:em_final],
      phi_shift_trace = phi_shift_hist[1:em_final],
      B_shift_trace = B_shift_hist[1:em_final],
      sigma2_shift_trace = sigma2_shift_hist[1:em_final],
      mu_r_hat_trace = mu_r_hat_trace[1:em_final],
      phi_r_hat_trace = phi_r_hat_trace[1:em_final],
      B_r_hat_trace = B_r_hat_trace[1:em_final],
      sigma2_r_hat_trace = sigma2_r_hat_trace[1:em_final],
      mu_remaining_trace = mu_remaining_trace[1:em_final],
      phi_remaining_trace = phi_remaining_trace[1:em_final],
      B_remaining_trace = B_remaining_trace[1:em_final],
      sigma2_remaining_trace = sigma2_remaining_trace[1:em_final],
      no_stall_trace = no_stall_trace[1:em_final],
      ev_ok_trace = ev_ok_trace[1:em_final]
    ))
  }

  out
}
