fit_pfa_woodbury_lam_corr <- function(Y, X, group, K,
                                      M = rowSums(Y),
                                      max_iter = 60, tol = 1e-4,
                                      lambda_phi = 0, sigma2_init = 0.3,
                                      verbose = FALSE,
                                      estep_max_iter = 100, estep_gtol = 1e-3,
                                      exact_Shat = FALSE,
                                      B_true = NULL,
                                      use_lambda_correction = FALSE,
                                      corr_max_rel = 0.5, 
                                      fix_sigma2 = NULL,    # NULL | numeric | "auto"
                                      freeze_tol = 1e-5, freeze_patience = 3,
                                      use_aitken = FALSE, aitken_window = 3, aitken_tol = 1e-4) {

  N <- nrow(Y); Q <- ncol(Y); P <- ncol(X); J <- max(group)
  stopifnot(all(X[, 1] == 1))
  stopifnot(is.null(fix_sigma2) || identical(fix_sigma2, "auto") || is.numeric(fix_sigma2))

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
  sigma2_trace <- numeric(max_iter)
  B_dist_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
  corr_rel_size_trace <- numeric(max_iter)
  corr_n_capped_trace <- integer(max_iter)
  Shat_share_trace <- numeric(max_iter)
  converged <- FALSE
  em <- 0

  sigma2_is_frozen <- FALSE; sigma2_freeze_iter <- NA_integer_; sigma2_freeze_value <- NA_real_
  freeze_stable_count <- 0
  if (is.numeric(fix_sigma2)) {
    sigma2 <- fix_sigma2
    sigma2_is_frozen <- TRUE; sigma2_freeze_iter <- 0L; sigma2_freeze_value <- fix_sigma2
    if (verbose) message(sprintf("sigma2 frozen from iter 0 at %.4f", fix_sigma2))
  }
  aitken_extrapolate <- function(x3) {
    d <- x3[3] - 2 * x3[2] + x3[1]
    if (!is.finite(d) || abs(d) < 1e-12) return(NA_real_)
    x3[3] - (x3[3] - x3[2])^2 / d
  }

  for (em in 1:max_iter) {

    es <- estep_wb(J, group, Y, X, M, mu, phi, B, sigma2,
                   max_iter = estep_max_iter, gtol = estep_gtol, exact_Shat = exact_Shat)
    lambda_hat <- es$lambda_hat
    S_hat <- es$S_hat

    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total
    
    n_capped <- 0L
    if (use_lambda_correction) {
      cc <- correct_lambda_edgeworth(lambda_hat, S_hat, Y, X, group, M, mu, phi)
      lam_c <- cc$lambda_corrected
      for (j in seq_len(J)) {                     # safety cap per group
        if (is.finite(cc$rel_size[j]) && cc$rel_size[j] > corr_max_rel) {
          mu1 <- lam_c[, j] - lambda_hat[, j]
          lam_c[, j] <- lambda_hat[, j] + mu1 * (corr_max_rel / cc$rel_size[j])
          n_capped <- n_capped + 1L
        }
      }
      lambda_hat <- lam_c
      corr_rel_size_trace[em] <- mean(pmin(cc$rel_size, corr_max_rel), na.rm = TRUE)
    } else {
      corr_rel_size_trace[em] <- 0
    }
    corr_n_capped_trace[em] <- n_capped

    mp <- mstep_phi_wb(Y, X, group, lambda_hat, mu, phi, lambda_phi = lambda_phi)
    mu <- mp$mu; phi <- mp$phi

    S_lam <- tcrossprod(lambda_hat) / J
    S_S <- matrix(0, Q, Q)
    for (j in 1:J) S_S <- S_S + S_hat[[j]] / J
    S_obs <- S_lam + S_S
    Shat_share_trace[em] <- sum(diag(S_S)) / sum(diag(S_obs))

    if (sigma2_is_frozen) {
      Ssym <- (S_obs + t(S_obs)) / 2
      eS <- eigen(Ssym, symmetric = TRUE)
      lamK <- eS$values[1:K]
      Uk <- eS$vectors[, 1:K, drop = FALSE]
      if (any(lamK < sigma2) && verbose)
        message(sprintf("iter %d: %d factor eigenvalue(s) below frozen sigma2 -- clipped", em, sum(lamK < sigma2)))
      B <- apply_PLT(Uk %*% diag(sqrt(pmax(lamK - sigma2, 0)), K))
    } else {
      rt <- rubin_thayer_wb(S_obs, K, B_init = B, sigma2_init = sigma2)
      B <- apply_PLT(rt$B)
      sigma2 <- rt$sigma2
    }

    sigma2_trace[em] <- sigma2
    if (!is.null(B_true)) B_dist_trace[em] <- subspace_dist(B, B_true)

    if (verbose) {
      bd <- if (!is.null(B_true)) sprintf("  B_dist=%.4f", B_dist_trace[em]) else ""
      cs <- if (use_lambda_correction) sprintf("  corr_sz=%.4f capped=%d", corr_rel_size_trace[em], n_capped) else ""
      fz <- if (sigma2_is_frozen) sprintf("  [FROZEN@%.4f]", sigma2_freeze_value) else ""
      cat(sprintf("iter %3d  log_ev=%.4f  sigma2=%.4f  Shat_share=%.3f%s%s%s\n",
                  em, log_ev[em], sigma2, Shat_share_trace[em], bd, cs, fz))
    }

    if (identical(fix_sigma2, "auto") && !sigma2_is_frozen) {
      if (use_aitken && em >= aitken_window + 1) {
        est_new <- aitken_extrapolate(sigma2_trace[(em - aitken_window + 1):em])
        est_old <- aitken_extrapolate(sigma2_trace[(em - aitken_window):(em - 1)])
        if (is.finite(est_new) && is.finite(est_old) &&
            abs(est_new - est_old) < aitken_tol * (abs(est_old) + 1e-8)) {
          sigma2 <- est_new
          sigma2_is_frozen <- TRUE; sigma2_freeze_iter <- em; sigma2_freeze_value <- est_new
          if (verbose) message(sprintf("iter %d: Aitken freeze at %.4f", em, est_new))
        }
      } else if (!use_aitken && em >= 2) {
        if (abs(sigma2_trace[em] - sigma2_trace[em - 1]) < freeze_tol) {
          freeze_stable_count <- freeze_stable_count + 1
        } else freeze_stable_count <- 0
        if (freeze_stable_count >= freeze_patience) {
          sigma2_is_frozen <- TRUE; sigma2_freeze_iter <- em; sigma2_freeze_value <- sigma2_trace[em]
          if (verbose) message(sprintf("iter %d: tol freeze at %.4f", em, sigma2_freeze_value))
        }
      }
    }

    if (em > 1 &&
        abs(log_ev[em] - log_ev[em - 1]) < tol * (abs(log_ev[em - 1]) + 1)) {
      converged <- TRUE; break
    }
  }

  list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
       lambda_hat = lambda_hat, S_hat = S_hat,
       log_evidence = log_ev[1:em], converged = converged, iterations = em,
       sigma2_trace = sigma2_trace[1:em],
       B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em] else NULL,
       corr_rel_size_trace = corr_rel_size_trace[1:em],
       corr_n_capped_trace = corr_n_capped_trace[1:em],
       Shat_share_trace = Shat_share_trace[1:em],
       sigma2_freeze_iter = sigma2_freeze_iter,
       sigma2_freeze_value = sigma2_freeze_value)
}
