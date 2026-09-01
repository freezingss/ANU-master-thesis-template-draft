library(SQUAREM)

pack_theta <- function(mu, phi, B, sigma2) {
  c(mu, as.vector(phi), as.vector(B), log(sigma2))
}

unpack_theta <- function(par, Q, P, K) {
  i1 <- Q
  i2 <- i1 + P * Q
  i3 <- i2 + Q * K
  list(mu = par[1:i1],
       phi = matrix(par[(i1 + 1):i2], P, Q),
       B = matrix(par[(i2 + 1):i3], Q, K),
       sigma2 = exp(par[i3 + 1]))
}

em_map_decorr <- function(par, Y, X, group, M, J, Q, P, K,
                          lambda_phi, estep_max_iter, estep_gtol,
                          use_lambda_correction, decouple_muphi_lambda,
                          corr_max_rel, env) {
  th <- unpack_theta(par, Q, P, K)
  mu <- th$mu
  phi <- th$phi
  B <- th$B
  sigma2 <- max(th$sigma2, 1e-8)

  es <- estep_wbonly(J, group, Y, X, M, mu, phi, B, sigma2,
                     estep_max_iter, estep_gtol)
  lambda_mode <- es$lambda_hat
  S_hat <- es$S_hat

  env$log_ev <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total

  if (use_lambda_correction) {
    cc <- correct_lambda_edgeworth(Q, J, lambda_mode, S_hat, Y, X, group, M, mu, phi)
    lam_c <- cc$lambda_corrected
    n_capped <- 0L
    for (j in seq_len(J)) {
      if (is.finite(cc$rel_size[j]) && cc$rel_size[j] > corr_max_rel) {
        d <- lam_c[, j] - lambda_mode[, j]
        lam_c[, j] <- lambda_mode[, j] + d * (corr_max_rel / cc$rel_size[j])
        n_capped <- n_capped + 1L
      }
    }
    lambda_corr <- lam_c
    env$rel_size <- mean(pmin(cc$rel_size, corr_max_rel), na.rm = TRUE)
    env$n_capped <- n_capped
  } else {
    lambda_corr <- lambda_mode
    env$rel_size <- 0
    env$n_capped <- 0L
  }

  lambda_for_muphi <- if (decouple_muphi_lambda) lambda_mode else lambda_corr

  mp <- mstep_phi_wbonly(Y, X, group, lambda_for_muphi, mu, phi, lambda_phi = lambda_phi)
  mu_new <- mp$mu
  phi_new <- mp$phi

  S_lam <- tcrossprod(lambda_corr) / J
  S_S <- matrix(0, Q, Q)
  for (j in 1:J) S_S <- S_S + S_hat[[j]] / J
  S_obs <- S_lam + S_S
  env$Shat_share <- sum(diag(S_S)) / sum(diag(S_obs))

  rt <- rubin_thayer_wbonly(S_obs, K, B_init = B, sigma2_init = sigma2)
  B_new <- apply_PLT(rt$B)
  sigma2_new <- max(rt$sigma2, 1e-8)

  env$lambda_mode <- lambda_mode
  env$lambda_corr <- lambda_corr
  env$S_hat <- S_hat

  pack_theta(mu_new, phi_new, B_new, sigma2_new)
}

neg_logev_decorr <- function(par, Y, X, group, M, J, Q, P, K,
                             lambda_phi, estep_max_iter, estep_gtol,
                             use_lambda_correction, decouple_muphi_lambda,
                             corr_max_rel, env) {
  th <- unpack_theta(par, Q, P, K)
  sigma2 <- max(th$sigma2, 1e-8)
  es <- tryCatch(estep_wbonly(J, group, Y, X, M, th$mu, th$phi, th$B, sigma2,
                              estep_max_iter, estep_gtol),
                 error = function(e) NULL)
  if (is.null(es)) return(Inf)
  val <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total
  if (!is.finite(val)) return(Inf)
  -val
}

fit_pfa_decorr_squarem <- function(Y, X, group, K,
                                   M = rowSums(Y),
                                   max_iter = 500, tol = 1e-3,
                                   lambda_phi = 0, sigma2_init = 0.3,
                                   verbose = FALSE,
                                   estep_max_iter = 100, estep_gtol = 1e-3,
                                   B_true = NULL,
                                   use_lambda_correction = TRUE,
                                   decouple_muphi_lambda = TRUE,
                                   corr_max_rel = 0.5,
                                   use_squarem = TRUE,
                                   objfn_inc = 0) {
  N <- nrow(Y)
  Q <- ncol(Y)
  P <- ncol(X)
  J <- max(group)
  stopifnot(all(X[, 1] == 1))

  avg_prop <- colMeans(Y / pmax(rowSums(Y), 1))
  mu0 <- log(avg_prop + 1e-8)
  mu0 <- mu0 - mean(mu0)
  phi0 <- matrix(0, P, Q)
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
  B0 <- apply_PLT(sv0$u %*% diag(pmax(sv0$d[1:K] * 0.5, 0.1), K))
  par0 <- pack_theta(mu0, phi0, B0, sigma2_init)

  env <- new.env(parent = emptyenv())
  env$log_ev <- NA_real_
  args <- list(Y = Y, X = X, group = group, M = M, J = J, Q = Q, P = P, K = K,
               lambda_phi = lambda_phi,
               estep_max_iter = estep_max_iter, estep_gtol = estep_gtol,
               use_lambda_correction = use_lambda_correction,
               decouple_muphi_lambda = decouple_muphi_lambda,
               corr_max_rel = corr_max_rel, env = env)

  t0 <- proc.time()["elapsed"]

  if (use_squarem) {
    ctrl <- list(tol = tol, maxiter = max_iter, objfn.inc = objfn_inc,
                 trace = verbose, intermed = FALSE)
    res <- do.call(squarem, c(list(par = par0,
                                   fixptfn = em_map_decorr,
                                   objfn = neg_logev_decorr,
                                   control = ctrl), args))
    par_final <- res$par
    n_fpeval <- res$fpevals
    n_objeval <- res$objfevals
    converged <- isTRUE(res$convergence)
    iters <- res$iter
    log_ev_trace <- NULL
  } else {
    par_cur <- par0
    log_ev_trace <- numeric(max_iter)
    converged <- FALSE
    iters <- 0
    for (t in 1:max_iter) {
      par_new <- do.call(em_map_decorr, c(list(par = par_cur), args))
      log_ev_trace[t] <- env$log_ev
      iters <- t
      if (t > 1 && abs(log_ev_trace[t] - log_ev_trace[t - 1]) < tol) {
        par_cur <- par_new
        converged <- TRUE
        break
      }
      par_cur <- par_new
    }
    par_final <- par_cur
    log_ev_trace <- log_ev_trace[1:iters]
    n_fpeval <- iters
    n_objeval <- 0
  }

  elapsed <- as.numeric(proc.time()["elapsed"] - t0)
  th <- unpack_theta(par_final, Q, P, K)
  final_logev <- -do.call(neg_logev_decorr, c(list(par = par_final), args))

  list(mu = th$mu, phi = th$phi, B = th$B, sigma2 = th$sigma2,
       lambda_mode = env$lambda_mode, lambda_corrected = env$lambda_corr,
       lambda_hat = env$lambda_corr, S_hat = env$S_hat,
       log_evidence_final = final_logev,
       log_evidence = log_ev_trace,
       converged = converged, iterations = iters,
       fpevals = n_fpeval, objfevals = n_objeval,
       time = elapsed,
       B_dist = if (!is.null(B_true)) subspace_dist(th$B, B_true) else NA_real_,
       Shat_share = env$Shat_share,
       accelerated = use_squarem)
}
