source("std_multi_wbonly.R")
source("existing_packages.R")

fit_pfa_wbonly_init_test <- function(Y, X, group, K,
                                     M = rowSums(Y),
                                     init_mode = c("warm", "cold", "warm_jitter"),
                                     init_jitter_sd = 0.1,
                                     max_iter = 500, tol = 1e-4,
                                     tol_mu = NULL, tol_phi = NULL,
                                     tol_sigma2 = NULL, tol_B = NULL,
                                     rel_tol = 1e-3,
                                     lambda_phi = 0, sigma2_init = 0.3,
                                     verbose = FALSE,
                                     estep_max_iter = 100, estep_gtol = 1e-3,
                                     trace = TRUE, B_true = NULL, sigma2_true = NULL,
                                     patience = 500, burn_in = 15, rate_window = 10,
                                     stall_patience = 500, also_require_mu_phi = FALSE,
                                     sigma2_floor = 1e-8) {
  
  init_mode <- match.arg(init_mode)
  N <- nrow(Y)
  Q <- ncol(Y)
  P <- ncol(X)
  J <- max(group)
  stopifnot(all(X[, 1] == 1))
  
  start_time <- Sys.time()
  
  if (init_mode == "cold") {
    mu <- rep(0, Q)
    phi <- matrix(0, P, Q)
    B <- matrix(0, Q, K)
  } else {
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
    
    if (init_mode == "warm_jitter") {
      mu <- mu + rnorm(Q, 0, init_jitter_sd)
      mu <- mu - mean(mu)
      B <- B + matrix(rnorm(Q * K, 0, init_jitter_sd), Q, K)
      B <- apply_PLT(B)
    }
  }
  sigma2 <- sigma2_init
  
  dt0 <- default_tols(mu, sigma2_init, K, rel_tol = rel_tol)
  if (is.null(tol_mu)) tol_mu <- dt0$tol_mu
  if (is.null(tol_phi)) tol_phi <- dt0$tol_phi
  if (is.null(tol_sigma2)) tol_sigma2 <- dt0$tol_sigma2
  if (is.null(tol_B)) tol_B <- dt0$tol_B
  
  lambda_hat <- matrix(0, Q, J)
  
  mu_init <- mu
  phi_init <- phi
  B_init <- B
  sigma2_init_val <- sigma2
  lambda_init_val <- lambda_hat
  
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
    tucker_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
    ang_max_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
    rv_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
    rmse_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
    sig_err_trace <- if (!is.null(B_true) && !is.null(sigma2_true)) numeric(max_iter) else NULL
    u_frac_trace <- numeric(max_iter)
    r_hat_trace <- rep(NA_real_, max_iter)
    remaining_gap_trace <- rep(NA_real_, max_iter)
    inner_iter_mean_trace <- rep(NA_real_, max_iter)
    grad_norm_mean_trace <- rep(NA_real_, max_iter)
    inner_maxiter_hit_pct_trace <- rep(NA_real_, max_iter)
    n_warnings_trace <- integer(max_iter)
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
    log_ev_delta_abs_trace <- rep(NA_real_, max_iter)
    log_ev_delta_rel_trace <- rep(NA_real_, max_iter)
    signal_eig_min_trace <- rep(NA_real_, max_iter)
    heywood_margin_trace <- rep(NA_real_, max_iter)
    sigma2_at_floor_trace <- rep(NA, max_iter)
    Shat_share_trace <- rep(NA_real_, max_iter)
    health_ok_trace <- rep(NA, max_iter)
    iter_wall_time_trace <- rep(NA_real_, max_iter)
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
    
    iter_start <- Sys.time()
    
    warn_msgs <- character(0)
    es <- withCallingHandlers(
      estep_wbonly_warm_start(J, group, Y, X, M, mu, phi, B, sigma2,
                              lambda_init = lambda_hat,
                              max_iter = estep_max_iter, estep_gtol = estep_gtol),
      warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)) })
    lambda_hat <- es$lambda_hat
    S_hat <- es$S_hat
    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total
    
    if (trace) monotone_trace[em] <- if (em == 1) TRUE else (log_ev[em] >= log_ev[em - 1])
    
    if (trace) {
      log_ev_delta_abs_trace[em] <- if (em == 1) NA_real_ else abs(log_ev[em] - log_ev[em - 1])
      log_ev_delta_rel_trace[em] <- if (em == 1) NA_real_ else log_ev_delta_abs_trace[em] / max(abs(log_ev[em - 1]), 1e-8)
      n_warnings_trace[em] <- length(warn_msgs)
    }
    
    if (log_ev[em] > best_log_ev) {
      best_log_ev <- log_ev[em]
      best_iter <- em
      best_state <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
                         lambda_hat = lambda_hat, S_hat = S_hat,
                         n_iter_vec = es$n_iter_vec, grad_norm_vec = es$grad_norm_vec)
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
    
    S_lam <- tcrossprod(lambda_hat) / J
    S_S <- matrix(0, Q, Q)
    for (j in 1:J) S_S <- S_S + S_hat[[j]] / J
    S_obs <- S_lam + S_S
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
      if (!is.null(B_true)) {
        B_dist_trace[em] <- subspace_distance(B, B_true)
        tucker_trace[em] <- tucker_congruence(B, B_true)$mean
        ang_max_trace[em] <- max(principal_angles(B, B_true))
        rv_trace[em] <- rv_coefficient(B, B_true)
        rmse_trace[em] <- loading_rmse(B, B_true)
        if (!is.null(sigma2_true)) sig_err_trace[em] <- sigma_error(B, sigma2, B_true, sigma2_true)
      }
      u_frac_trace[em] <- u_frac(B)
      inner_iter_mean_trace[em] <- mean(es$n_iter_vec, na.rm = TRUE)
      grad_norm_mean_trace[em] <- mean(es$grad_norm_vec, na.rm = TRUE)
      inner_maxiter_hit_pct_trace[em] <- mean(es$n_iter_vec >= estep_max_iter, na.rm = TRUE)
      hd <- heywood_diagnostics(rt$ev, K, sigma2)
      signal_eig_min_trace[em] <- hd$signal_eig_min
      heywood_margin_trace[em] <- hd$heywood_margin
      sigma2_at_floor_trace[em] <- sigma2 <= sigma2_floor + 1e-10
      Shat_share_trace[em] <- sum(diag(S_S)) / sum(diag(S_obs))
      health_ok_trace[em] <- is.finite(log_ev[em]) && is.finite(sigma2) && all(is.finite(B)) && all(is.finite(mu))
      iter_wall_time_trace[em] <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
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
  
  end_time <- Sys.time()
  
  final_literal <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
                        lambda_hat = lambda_hat, S_hat = S_hat,
                        ev = rt$ev, n_iter_vec = es$n_iter_vec, grad_norm_vec = es$grad_norm_vec)
  
  out <- list(mu = best_state$mu, phi = best_state$phi,
              B = best_state$B, sigma2 = best_state$sigma2,
              lambda_hat = best_state$lambda_hat, S_hat = best_state$S_hat,
              log_evidence = log_ev[1:em], best_iter = best_iter, best_log_ev = best_log_ev,
              converged = converged, iterations = em, stop_reason = stop_reason,
              init_mode = init_mode,
              tol = tol, tol_mu = tol_mu, tol_phi = tol_phi, tol_sigma2 = tol_sigma2, tol_B = tol_B,
              sigma2_floor = sigma2_floor,
              meta = list(pbs_job_id = Sys.getenv("PBS_JOBID"),
                          start_time = start_time, end_time = end_time,
                          elapsed_sec = as.numeric(difftime(end_time, start_time, units = "secs"))),
              checkpoints = list(
                init = list(mu = mu_init, phi = phi_init, B = B_init,
                            sigma2 = sigma2_init_val, lambda_hat = lambda_init_val),
                peak = best_state,
                final = final_literal))
  
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
      tucker_trace = if (!is.null(B_true)) tucker_trace[1:em_final] else NULL,
      ang_max_trace = if (!is.null(B_true)) ang_max_trace[1:em_final] else NULL,
      rv_trace = if (!is.null(B_true)) rv_trace[1:em_final] else NULL,
      rmse_trace = if (!is.null(B_true)) rmse_trace[1:em_final] else NULL,
      sig_err_trace = if (!is.null(B_true) && !is.null(sigma2_true)) sig_err_trace[1:em_final] else NULL,
      u_frac_trace = u_frac_trace[1:em_final],
      r_hat_trace = r_hat_trace[1:em_final],
      remaining_gap_trace = remaining_gap_trace[1:em_final],
      inner_iter_mean_trace = inner_iter_mean_trace[1:em_final],
      grad_norm_mean_trace = grad_norm_mean_trace[1:em_final],
      inner_maxiter_hit_pct_trace = inner_maxiter_hit_pct_trace[1:em_final],
      n_warnings_trace = n_warnings_trace[1:em_final],
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
      ev_ok_trace = ev_ok_trace[1:em_final],
      log_ev_delta_abs_trace = log_ev_delta_abs_trace[1:em_final],
      log_ev_delta_rel_trace = log_ev_delta_rel_trace[1:em_final],
      signal_eig_min_trace = signal_eig_min_trace[1:em_final],
      heywood_margin_trace = heywood_margin_trace[1:em_final],
      sigma2_at_floor_trace = sigma2_at_floor_trace[1:em_final],
      Shat_share_trace = Shat_share_trace[1:em_final],
      health_ok_trace = health_ok_trace[1:em_final],
      iter_wall_time_trace = iter_wall_time_trace[1:em_final]
    ))
  }
  
  out
}

fit_wb_cold <- function(Y, X, group, K, ...) {
  fit_pfa_wbonly_init_test(Y, X, group, K, init_mode = "cold", ...)
}

fit_wb_warm <- function(Y, X, group, K, ...) {
  fit_pfa_wbonly_init_test(Y, X, group, K, init_mode = "warm", ...)
}

fit_wb_warm_jitter <- function(Y, X, group, K, init_jitter_sd = 0.1, ...) {
  fit_pfa_wbonly_init_test(Y, X, group, K, init_mode = "warm_jitter",
                           init_jitter_sd = init_jitter_sd, ...)
}

fit_glmm_multistart <- function(Y, X, group, K, n_restart = 5, jitter_sd = 0.2) {
  overall_start <- Sys.time()
  long <- build_long(Y, X, group)
  xnames <- paste0("x", 2:ncol(X))
  covar_terms <- paste(paste0("category:", xnames), collapse = " + ")
  form <- as.formula(paste0("count ~ category + ", covar_terms,
                            " + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))
  ctrl <- glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  
  extract <- function(fit, elapsed_sec) {
    conv <- isTRUE(fit$sdr$pdHess)
    L <- tryCatch({
      vc <- as.matrix(glmmTMB::VarCorr(fit)$cond$group)
      e <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
      e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(e$values[1:K], 0)), K)
    }, error = function(e) NULL)
    if (is.null(L)) {
      return(list(B = NULL, sigma2 = NA_real_, ok = FALSE, converged = conv,
                  loglik = NA_real_, aic = NA_real_, error = "loading extraction failed",
                  elapsed_sec = elapsed_sec))
    }
    sigma2_hat <- tryCatch(as.numeric(glmmTMB::VarCorr(fit)$cond$obs)[1], error = function(e) NA_real_)
    loglik_hat <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
    aic_hat <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
    list(B = L, sigma2 = sigma2_hat, ok = TRUE, converged = conv,
         loglik = loglik_hat, aic = aic_hat, error = NA_character_,
         elapsed_sec = elapsed_sec)
  }
  fit0_start <- Sys.time()
  fit0 <- tryCatch(glmmTMB::glmmTMB(form, data = long, family = poisson(link = "log"),
                                    offset = log_total, control = ctrl),
                   error = function(e) e)
  fit0_elapsed <- as.numeric(difftime(Sys.time(), fit0_start, units = "secs"))
  if (inherits(fit0, "error")) {
    fail <- list(B = NULL, sigma2 = NA_real_, ok = FALSE, converged = FALSE,
                 loglik = NA_real_, aic = NA_real_, error = conditionMessage(fit0),
                 elapsed_sec = fit0_elapsed)
    overall_elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
    return(list(best = fail, best_idx = NA_integer_, all_restarts = rep(list(fail), n_restart),
                meta = list(pbs_job_id = Sys.getenv("PBS_JOBID"), elapsed_sec = overall_elapsed)))
  }
  
  restarts <- vector("list", n_restart)
  restarts[[1]] <- extract(fit0, fit0_elapsed)
  
  if (n_restart > 1) {
    pl0 <- fit0$obj$env$parList()
    rm(fit0)
    gc(verbose = FALSE)
    for (r in 2:n_restart) {
      r_start <- Sys.time()
      pl_jit <- pl0
      pl_jit$b <- pl0$b + rnorm(length(pl0$b), 0, jitter_sd)
      pl_jit$beta <- pl0$beta + rnorm(length(pl0$beta), 0, jitter_sd)
      fit_r <- tryCatch(glmmTMB::glmmTMB(form, data = long, family = poisson(link = "log"),
                                         offset = log_total, control = ctrl, start = pl_jit),
                        error = function(e) e)
      r_elapsed <- as.numeric(difftime(Sys.time(), r_start, units = "secs"))
      restarts[[r]] <- if (inherits(fit_r, "error")) {
        list(B = NULL, sigma2 = NA_real_, ok = FALSE, converged = FALSE,
             loglik = NA_real_, aic = NA_real_, error = conditionMessage(fit_r),
             elapsed_sec = r_elapsed)
      } else extract(fit_r, r_elapsed)
      rm(fit_r)
      gc(verbose = FALSE)
    }
  }
  
  logliks <- sapply(restarts, function(r) if (isTRUE(r$ok) && is.finite(r$loglik)) r$loglik else -Inf)
  best_idx <- if (all(!is.finite(logliks))) NA_integer_ else which.max(logliks)
  best <- if (is.na(best_idx)) restarts[[1]] else restarts[[best_idx]]
  
  overall_elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
  
  list(best = best, best_idx = best_idx, all_restarts = restarts,
       meta = list(pbs_job_id = Sys.getenv("PBS_JOBID"), elapsed_sec = overall_elapsed))
}