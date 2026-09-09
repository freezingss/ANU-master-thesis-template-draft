source("basic_functions.R")
source("pfa_woodbury_only.R")

deconv_sigma <- function(es, tau2, Q, J, ridge = 1e-3) {
  Sig <- matrix(0, Q, Q)
  nbad <- 0L
  for (j in 1:J) {
    Sj <- es$S_hat[[j]]
    Hj <- tryCatch(chol2inv(chol((Sj + t(Sj)) / 2)), error = function(e) NULL)
    if (is.null(Hj)) { nbad <- nbad + 1L; next }
    Fj <- Hj - diag(Q) / tau2
    Fr <- (Fj + t(Fj)) / 2 + ridge * diag(Q)
    Fi <- tryCatch(chol2inv(chol(Fr)), error = function(e) NULL)
    if (is.null(Fi)) { nbad <- nbad + 1L; next }
    lt <- es$lambda_hat[, j] + as.numeric(Fi %*% es$lambda_hat[, j]) / tau2
    Sig <- Sig + (tcrossprod(lt) - Fi) / J
  }
  if (nbad > 0) warning(sprintf("deconv_sigma: %d of %d groups failed (Cholesky) and were skipped", nbad, J))
  list(S = (Sig + t(Sig)) / 2, n_failed = nbad)
}

fit_pfa_debias <- function(Y, X, group, K,
                           M = rowSums(Y),
                           max_iter = 100, tol = 1e-4,
                           lambda_phi = 0, sigma2_init = 0.3,
                           estep_max_iter = 100, estep_gtol = 1e-3,
                           ridge = 1e-3, verbose = FALSE, B_true = NULL,
                           trace = TRUE) {
  Q <- ncol(Y); J <- max(group)
  stopifnot(all(X[, 1] == 1))
  th <- init_theta(Y, X, group, K, sigma2_init)
  mu <- th$mu; phi <- th$phi; B <- th$B; sigma2 <- th$sigma2
  tau2 <- sum(diag(tcrossprod(B))) / Q + sigma2
  
  log_ev <- numeric(max_iter); s2_tr <- numeric(max_iter)
  monotone_trace <- rep(NA, max_iter)
  B_dist_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
  converged <- FALSE; em <- 0; nbad <- 0L
  
  if (trace) {
    best_iter <- NA_integer_
    best_log_ev <- -Inf
    best_state <- NULL
  }
  
  for (em in 1:max_iter) {
    es <- estep_wbonly(J, group, Y, X, M, mu, phi, matrix(0, Q, K), tau2,
                       estep_max_iter, estep_gtol)
    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total
    monotone_trace[em] <- if (em == 1) TRUE else (log_ev[em] >= log_ev[em - 1])
    
    if (trace && log_ev[em] > best_log_ev) {
      best_log_ev <- log_ev[em]
      best_iter <- em
      best_state <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2, tau2 = tau2,
                         lambda_hat = es$lambda_hat, S_hat = es$S_hat)
    }
    
    mp <- mstep_phi_wbonly(Y, X, group, es$lambda_hat, mu, phi, lambda_phi = lambda_phi)
    mu <- mp$mu; phi <- mp$phi
    
    dc <- deconv_sigma(es, tau2, Q, J, ridge)
    nbad <- nbad + dc$n_failed
    pc <- ppca_closed(dc$S, K)
    B <- apply_PLT(pc$B); sigma2 <- pc$sigma2
    tau2 <- sum(diag(tcrossprod(B))) / Q + sigma2
    s2_tr[em] <- sigma2
    if (!is.null(B_true)) B_dist_trace[em] <- subspace_distance(B, B_true)
    
    if (verbose) {
      bd <- if (!is.null(B_true)) sprintf(" B_dist=%.4f", B_dist_trace[em]) else ""
      cat(sprintf("iter %3d  log_ev=%.4f  sigma2=%.4f  tau2=%.4f%s\n",
                  em, log_ev[em], sigma2, tau2, bd))
    }
    if (em > 1 && abs(log_ev[em] - log_ev[em - 1]) < tol) { converged <- TRUE; break }
  }
  
  es_final <- estep_wbonly(J, group, Y, X, M, mu, phi, matrix(0, Q, K), tau2,
                           estep_max_iter, estep_gtol)
  log_ev_final <- es_final$lp_total - 0.5 * J * es_final$log_det_Sigma + 0.5 * es_final$ld_S_total
  
  em_final <- em
  is_monotone <- all(monotone_trace[1:em_final])
  first_decrease_iter <- if (!is_monotone) which(!monotone_trace[1:em_final])[1] else NA_integer_
  
  out <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
              lambda_hat = es_final$lambda_hat, S_hat = es_final$S_hat,
              log_evidence = log_ev[1:em_final], log_ev_final = log_ev_final,
              sigma2_trace = s2_tr[1:em_final],
              monotone_trace = monotone_trace[1:em_final], is_monotone = is_monotone,
              first_decrease_iter = first_decrease_iter,
              B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em_final] else NULL,
              converged = converged, iterations = em_final, n_failed_groups = nbad, tau2 = tau2)
  
  if (trace) {
    out <- c(out, list(best_iter = best_iter, best_log_ev = best_log_ev, best = best_state))
  }
  out
}

fit_pfa_decorr <- function(Y, X, group, K,
                           M = rowSums(Y),
                           max_iter = 100, tol = 1e-4,
                           lambda_phi = 0, sigma2_init = 0.3,
                           estep_max_iter = 100, estep_gtol = 1e-3,
                           corr_max_rel = 0.5, use_ppca = FALSE, verbose = FALSE,
                           B_true = NULL, trace = TRUE) {
  Q <- ncol(Y); J <- max(group)
  th <- init_theta(Y, X, group, K, sigma2_init)
  mu <- th$mu; phi <- th$phi; B <- th$B; sigma2 <- th$sigma2
  
  log_ev <- numeric(max_iter); converged <- FALSE; em <- 0
  sigma2_trace <- numeric(max_iter)
  monotone_trace <- rep(NA, max_iter)
  B_dist_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
  corr_rel_size_trace <- numeric(max_iter)
  corr_n_capped_trace <- integer(max_iter)
  
  if (trace) {
    best_iter <- NA_integer_
    best_log_ev <- -Inf
    best_state <- NULL
  }
  
  for (em in 1:max_iter) {
    es <- estep_wbonly(J, group, Y, X, M, mu, phi, B, sigma2, estep_max_iter, estep_gtol)
    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total
    monotone_trace[em] <- if (em == 1) TRUE else (log_ev[em] >= log_ev[em - 1])
    
    ac <- apply_lambda_correction(Q, J, es$lambda_hat, es$S_hat, Y, X, group, M, mu, phi, corr_max_rel)
    lc <- ac$lambda_corrected
    corr_rel_size_trace[em] <- mean(pmin(ac$rel_size, corr_max_rel), na.rm = TRUE)
    corr_n_capped_trace[em] <- ac$n_capped
    
    if (trace && log_ev[em] > best_log_ev) {
      best_log_ev <- log_ev[em]
      best_iter <- em
      best_state <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
                         lambda_hat = lc, S_hat = es$S_hat)
    }
    
    mp <- mstep_phi_wbonly(Y, X, group, lc, mu, phi, lambda_phi = lambda_phi)
    mu <- mp$mu; phi <- mp$phi
    
    SS <- matrix(0, Q, Q)
    for (j in 1:J) SS <- SS + es$S_hat[[j]] / J
    S_obs <- tcrossprod(lc) / J + SS
    if (use_ppca) {
      pc <- ppca_closed(S_obs, K); B <- apply_PLT(pc$B); sigma2 <- pc$sigma2
    } else {
      rt <- rubin_thayer_wbonly(S_obs, K, B_init = B, sigma2_init = sigma2)
      B <- apply_PLT(rt$B); sigma2 <- rt$sigma2
    }
    sigma2_trace[em] <- sigma2
    if (!is.null(B_true)) B_dist_trace[em] <- subspace_distance(B, B_true)
    
    if (verbose) {
      bd <- if (!is.null(B_true)) sprintf("  B_dist=%.4f", B_dist_trace[em]) else ""
      cat(sprintf("iter %3d  log_ev=%.4f  sigma2=%.4f%s\n", em, log_ev[em], sigma2, bd))
    }
    if (em > 1 && abs(log_ev[em] - log_ev[em - 1]) < tol) { converged <- TRUE; break }
  }
  
  es_final <- estep_wbonly(J, group, Y, X, M, mu, phi, B, sigma2, estep_max_iter, estep_gtol)
  ac_final <- apply_lambda_correction(Q, J, es_final$lambda_hat, es_final$S_hat, Y, X, group, M, mu, phi, corr_max_rel)
  lambda_hat <- ac_final$lambda_corrected
  S_hat <- es_final$S_hat
  log_ev_final <- es_final$lp_total - 0.5 * J * es_final$log_det_Sigma + 0.5 * es_final$ld_S_total
  
  em_final <- em
  is_monotone <- all(monotone_trace[1:em_final])
  first_decrease_iter <- if (!is_monotone) which(!monotone_trace[1:em_final])[1] else NA_integer_
  
  out <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
              lambda_hat = lambda_hat, S_hat = S_hat,
              log_evidence = log_ev[1:em_final], log_ev_final = log_ev_final,
              sigma2_trace = sigma2_trace[1:em_final],
              monotone_trace = monotone_trace[1:em_final], is_monotone = is_monotone,
              first_decrease_iter = first_decrease_iter,
              B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em_final] else NULL,
              corr_rel_size_trace = corr_rel_size_trace[1:em_final],
              corr_n_capped_trace = corr_n_capped_trace[1:em_final],
              converged = converged, iterations = em_final)
  
  if (trace) {
    out <- c(out, list(best_iter = best_iter, best_log_ev = best_log_ev, best = best_state))
  }
  out
}

fit_glmmTMB_ref <- function(Y, X, group, K) {
  long <- build_long(Y, X, group)
  form <- as.formula(paste0("count ~ category + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))
  ctrl <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  fit <- tryCatch(glmmTMB(form, data = long, family = poisson(link = "log"),
                          offset = long$log_total, control = ctrl),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    ctrl2 <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000),
                            start_method = list(method = "res"))
    fit <- tryCatch(glmmTMB(form, data = long, family = poisson(link = "log"),
                            offset = long$log_total, control = ctrl2),
                    error = function(e) e)
    if (inherits(fit, "error")) return(list(B = NULL, sigma2 = NA_real_, ok = FALSE,
                                            converged = FALSE))
  }
  L <- tryCatch({
    vc <- as.matrix(VarCorr(fit)$cond$group)
    e <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
    e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(e$values[1:K], 0)), K)
  }, error = function(e) NULL)
  if (is.null(L)) return(list(B = NULL, sigma2 = NA_real_, ok = FALSE,
                              converged = isTRUE(fit$sdr$pdHess)))
  list(B = L, sigma2 = tryCatch(as.numeric(VarCorr(fit)$cond$obs)[1],
                                error = function(e) NA_real_),
       ok = TRUE, converged = isTRUE(fit$sdr$pdHess))
}

# run_grid <- function(Js, Qs, seeds, K = 2, Nj = 15, sigma2_true = 0.3, P = 3,
#                      M_rate = 150, max_iter = 100, tol = 1e-4,
#                      methods = c("glmmTMB", "wb", "decorr", "debias"),
#                      verbose = TRUE) {
#   out <- NULL
#   for (J in Js) for (Q in Qs) for (sd in seeds) {
#     dat <- simulate_pfa_data(Q = Q, K = K, J = J, P = P, N_per_group = Nj,
#                              sigma2 = sigma2_true, M_rate = M_rate, seed = sd)
#     Bt <- dat$true$B
#     fl <- bbp_floor(Bt, sigma2_true, J)
#     orc <- ppca_closed(tcrossprod(dat$true$lambda) / J, K)
#     rows <- list()
#     rows[["oracle"]] <- list(t = 0, B = orc$B, s2 = orc$sigma2, cv = NA, it = 0)
#     for (m in methods) {
#       if (m == "glmmTMB") {
#         tt <- system.time(f <- fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K))[3]
#         rows[[m]] <- list(t = as.numeric(tt), B = f$B, s2 = f$sigma2,
#                           cv = f$converged, it = NA)
#       } else if (m == "wb") {
#         tt <- system.time(f <- fit_pfa_wbonly(dat$Y, dat$X, dat$group, K = K,
#                                               M = dat$M, max_iter = max_iter, tol = tol, sigma2_init = sigma2_true,
#                                               verbose = FALSE, estep_max_iter = 100, estep_gtol = 1e-3,
#                                               trace = FALSE))[3]
#         rows[[m]] <- list(t = as.numeric(tt), B = f$B, s2 = f$sigma2,
#                           cv = f$converged, it = f$iterations)
#       } else if (m == "decorr") {
#         tt <- system.time(f <- fit_pfa_decorr(dat$Y, dat$X, dat$group, K, M = dat$M,
#                                               max_iter = max_iter, tol = tol, sigma2_init = sigma2_true,
#                                               trace = FALSE))[3]
#         rows[[m]] <- list(t = as.numeric(tt), B = f$B, s2 = f$sigma2,
#                           cv = f$converged, it = f$iterations)
#       } else if (m == "debias") {
#         tt <- system.time(f <- fit_pfa_debias(dat$Y, dat$X, dat$group, K, M = dat$M,
#                                               max_iter = max_iter, tol = tol, sigma2_init = sigma2_true,
#                                               trace = FALSE))[3]
#         rows[[m]] <- list(t = as.numeric(tt), B = f$B, s2 = f$sigma2,
#                           cv = f$converged, it = f$iterations)
#       }
#     }
#     for (nm in names(rows)) {
#       r <- rows[[nm]]
#       mt <- B_metrics(r$B, r$s2, Bt, sigma2_true)
#       out <- rbind(out, data.frame(J = J, Q = Q, seed = sd, method = nm,
#                                    time_s = r$t, sigma2 = r$s2, d_true = mt["d_true"], ang_max = mt["ang_max"],
#                                    tucker = mt["tucker"], rmse = mt["rmse"], sig_err = mt["sig_err"],
#                                    floor = fl, converged = r$cv, iters = r$it,
#                                    stringsAsFactors = FALSE, row.names = NULL))
#     }
#     if (verbose) {
#       cat(sprintf("\n=== J=%d Q=%d seed=%d | true sigma2=%.3f | BBP floor=%.4f ===\n",
#                   J, Q, sd, sigma2_true, fl))
#       cat(sprintf("%-9s %8s %8s %8s %8s %8s %8s %8s %6s\n", "method", "time",
#                   "sigma2", "d_true", "ang_max", "tucker", "rmse", "sig_err", "conv"))
#       sub <- out[out$J == J & out$Q == Q & out$seed == sd, ]
#       for (i in seq_len(nrow(sub))) {
#         f4 <- function(x) if (is.na(x)) "-" else sprintf("%.4f", x)
#         cat(sprintf("%-9s %8.1f %8s %8s %8s %8s %8s %8s %6s\n",
#                     sub$method[i], sub$time_s[i], f4(sub$sigma2[i]), f4(sub$d_true[i]),
#                     if (is.na(sub$ang_max[i])) "-" else sprintf("%.1f", sub$ang_max[i]),
#                     f4(sub$tucker[i]), f4(sub$rmse[i]), f4(sub$sig_err[i]),
#                     if (is.na(sub$converged[i])) "-" else ifelse(sub$converged[i], "yes", "NO")))
#       }
#     }
#   }
#   out
# }
# 
# summarize_grid <- function(res) {
#   vs <- c("time_s", "sigma2", "d_true", "ang_max", "tucker", "rmse", "sig_err")
#   ag <- Reduce(function(a, b) merge(a, b, by = c("method", "J", "Q")),
#                lapply(vs, function(v)
#                  aggregate(as.formula(paste(v, "~ method + J + Q")), data = res,
#                            FUN = function(x) mean(x, na.rm = TRUE))))
#   ag <- ag[order(ag$J, ag$Q, ag$method), ]
#   cat("\n===== averaged over seeds =====\n")
#   cat(sprintf("%-9s %5s %5s | %7s %8s %8s %8s %8s %8s %8s\n", "method", "J", "Q",
#               "time", "sigma2", "d_true", "ang_max", "tucker", "rmse", "sig_err"))
#   for (i in seq_len(nrow(ag))) {
#     cat(sprintf("%-9s %5d %5d | %7.1f %8.4f %8.4f %8.1f %8.4f %8.4f %8.4f\n",
#                 ag$method[i], ag$J[i], ag$Q[i], ag$time_s[i], ag$sigma2[i],
#                 ag$d_true[i], ag$ang_max[i], ag$tucker[i], ag$rmse[i], ag$sig_err[i]))
#   }
#   pw <- NULL
#   for (J in unique(res$J)) for (Q in unique(res$Q)) {
#     a <- res[res$J == J & res$Q == Q & res$method == "decorr", ]
#     b <- res[res$J == J & res$Q == Q & res$method == "debias", ]
#     if (!nrow(a) || !nrow(b)) next
#     m <- merge(a, b, by = "seed", suffixes = c(".dec", ".deb"))
#     d <- m$d_true.deb - m$d_true.dec
#     pw <- rbind(pw, data.frame(J = J, Q = Q, n = length(d), mean_diff = mean(d),
#                                se = sd(d) / sqrt(length(d)), n_better = sum(d < 0)))
#   }
#   if (!is.null(pw)) {
#     cat("\n===== paired d_true difference (debias - decorr), negative favours debias =====\n")
#     cat(sprintf("%5s %5s %4s %11s %10s %10s\n", "J", "Q", "n", "mean_diff", "se", "n_better"))
#     for (i in seq_len(nrow(pw)))
#       cat(sprintf("%5d %5d %4d %11.4f %10.4f %10d\n", pw$J[i], pw$Q[i], pw$n[i],
#                   pw$mean_diff[i], pw$se[i], pw$n_better[i]))
#   }
#   invisible(list(avg = ag, paired = pw))
# }