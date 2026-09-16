source("basic_functions.R")
source("pfa_base_dense.R")
source("pfa_woodbury_only.R")
source("pfa_woodbury_lam_corrected.R")
source("pfa_debias.R")
source("pfa_fhem.R")
source("pfa_fhem_woodbury.R")
source("sim_data.R")

res_grid <- run_grid(Js = c(50, 75, 100),
                     Qs = c(75, 100, 125, 150),
                     seeds = 1,
                     K = 2, Nj = 15, sigma2_true = 0.3,
                     M_range = c(75, 150),
                     max_iter = 80, tol = 1e-4,
                     methods = c("wb", "corrected", "fhem", "glmmTMB"))
saveRDS(res_grid, file = "res_grid_multiJQ.rds")

res_reps <- run_grid(Js = c(10),
                     Qs = c(10),
                     seeds = 1,
                     K = 2, 
                     Nj = 15, 
                     sigma2_true = 0.3,
                     M_range = c(10, 20),
                     max_iter = 80, 
                     tol = 1e-4,
                     methods = c("wb", "corrected", "fhem", "glmmTMB"))
# saveRDS(res_reps, file = "res_reps_singleJQ4.rds")

plot_log_evidence(attr(res_reps, "traj_store"))
plot_B_trajectory(attr(res_reps, "traj_store"))
plot_sigma2_trajectory(attr(res_reps, "traj_store"), sigma2_true = 0.3)
plot_u_frac_trajectory(attr(res_reps, "traj_store"))
plot_runtime(res_reps)

plot_runtime(res_grid)

# Run grid function
run_grid <- function(Js, Qs, seeds, K = 2, Nj = 15, sigma2_true = 0.3, P = 3,
                     N_per_group_range = NULL, M_range = c(150, 150),
                     dispersion = Inf,
                     outlier_frac = 0, outlier_type = c("spike", "random_composition", "shock"),
                     outlier_severity = NULL,
                     max_iter = 100, tol = 1e-4,
                     methods = c("wb", "corrected", "fhem", "glmmTMB"),
                     verbose = TRUE) {
  outlier_type <- match.arg(outlier_type)
  if (is.null(N_per_group_range)) N_per_group_range <- c(Nj, Nj)
  out <- NULL
  B_store <- list()
  traj_store <- list()
  for (J in Js) for (Q in Qs) for (sd in seeds) {
    dat <- simulate_pfa_data(Q = Q, K = K, J = J, P = P,
                             N_per_group_range = N_per_group_range,
                             sigma2 = sigma2_true, M_range = M_range,
                             dispersion = dispersion,
                             outlier_frac = outlier_frac, outlier_type = outlier_type,
                             outlier_severity = outlier_severity,
                             seed = sd)
    Bt <- dat$true$B
    fl <- bbp_floor(Bt, sigma2_true, J)
    Pi_lam <- dat$true$lambda - matrix(colMeans(dat$true$lambda), nrow(dat$true$lambda),
                                       ncol(dat$true$lambda), byrow = TRUE)
    orc <- ppca_closed(tcrossprod(Pi_lam) / J, K)
    rows <- list()
    rows[["oracle"]] <- list(t = 0, B = orc$B, s2 = orc$sigma2, cv = NA, it = 0, traj = NULL)
    
    for (m in methods) {
      if (m == "glmmTMB") {
        tt <- system.time(f <- fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K))[3]
        rows[[m]] <- list(t = as.numeric(tt), B = f$B, s2 = f$sigma2,
                          cv = f$converged, it = NA, traj = NULL)
      } else if (m == "wb") {
        tt <- system.time(f <- fit_pfa_wbonly(dat$Y, dat$X, dat$group, K = K,
                                              M = dat$M, max_iter = max_iter, tol = tol, sigma2_init = sigma2_true,
                                              verbose = FALSE, estep_max_iter = 100, estep_gtol = 1e-3,
                                              trace = TRUE, B_true = Bt))[3]
        rows[[m]] <- list(t = as.numeric(tt), B = f$B, s2 = f$sigma2,
                          cv = f$converged, it = f$iterations, traj = f)
      } else if (m == "corrected") {
        tt <- system.time(f <- fit_pfa_woodbury_lam_corr(dat$Y, dat$X, dat$group, K,
                                                         M = dat$M, max_iter = max_iter, tol = tol, sigma2_init = sigma2_true,
                                                         verbose = FALSE, estep_max_iter = 100, estep_gtol = 1e-3,
                                                         use_lambda_correction = TRUE, corr_max_rel = 0.5,
                                                         trace = TRUE, B_true = Bt))[3]
        rows[[m]] <- list(t = as.numeric(tt), B = f$B, s2 = f$sigma2,
                          cv = f$converged, it = f$iterations, traj = f)
      } else if (m == "fhem") {
        tt <- system.time(f <- fit_pfa_fhem(dat$Y, dat$X, dat$group, K, M = dat$M,
                                            max_iter = 30, tol = tol, sigma2_init = sigma2_true,
                                            verbose = FALSE, estep_max_iter = 100, estep_gtol = 1e-3,
                                            fhem_max_iter = 300, fhem_tol = 1e-6,
                                            ridge = 0, deflate = TRUE, trace = TRUE, B_true = Bt))[3]
        rows[[m]] <- list(t = as.numeric(tt), B = f$B, s2 = f$sigma2,
                          cv = f$converged, it = f$iterations, traj = f)
      }
    }
    
    Bt_d <- deflate_cols(Bt)
    for (nm in names(rows)) {
      r <- rows[[nm]]
      mt_raw <- B_metrics(r$B, r$s2, Bt, sigma2_true)
      mt <- B_metrics(deflate_cols(r$B), r$s2, Bt_d, sigma2_true)
      out <- rbind(out, data.frame(J = J, Q = Q, seed = sd, method = nm,
                                   time_s = r$t, sigma2 = r$s2, u_frac = u_frac(r$B),
                                   d_true = mt["d_true"], ang_max = mt["ang_max"],
                                   tucker = mt["tucker"], rv = mt["rv"], rmse = mt["rmse"], sig_err = mt["sig_err"],
                                   d_true_raw = mt_raw["d_true"], tucker_raw = mt_raw["tucker"],
                                   rmse_raw = mt_raw["rmse"], sig_err_raw = mt_raw["sig_err"],
                                   floor = fl, converged = r$cv, iters = r$it,
                                   stringsAsFactors = FALSE, row.names = NULL))
      B_store[[length(B_store) + 1]] <- list(J = J, Q = Q, seed = sd, method = nm,
                                             B_hat = r$B, B_hat_deflated = deflate_cols(r$B),
                                             B_true = Bt, B_true_deflated = Bt_d)
      if (!is.null(r$traj)) {
        traj_store[[length(traj_store) + 1]] <- list(
          J = J, Q = Q, seed = sd, method = nm, floor = fl,
          log_evidence = r$traj$log_evidence,
          sigma2_trace = r$traj$sigma2_trace,
          B_dist_trace = r$traj$B_dist_trace,
          u_frac_trace = r$traj$u_frac_trace)
      }
    }
    if (verbose) {
      cat(sprintf("\n=== J=%d Q=%d seed=%d | true sigma2=%.3f | BBP floor=%.4f ===\n",
                  J, Q, sd, sigma2_true, fl))
      cat(sprintf("%-9s %8s %8s %8s %8s %8s %8s %8s %8s %6s\n", "method", "time",
                  "sigma2", "d_true", "ang_max", "tucker", "rv", "rmse", "sig_err", "conv"))
      sub <- out[out$J == J & out$Q == Q & out$seed == sd, ]
      for (i in seq_len(nrow(sub))) {
        f4 <- function(x) if (is.na(x)) "-" else sprintf("%.4f", x)
        cat(sprintf("%-9s %8.1f %8s %8s %8s %8s %8s %8s %8s %6s\n",
                    sub$method[i], sub$time_s[i], f4(sub$sigma2[i]), f4(sub$d_true[i]),
                    if (is.na(sub$ang_max[i])) "-" else sprintf("%.1f", sub$ang_max[i]),
                    f4(sub$tucker[i]), f4(sub$rv[i]), f4(sub$rmse[i]), f4(sub$sig_err[i]),
                    if (is.na(sub$converged[i])) "-" else ifelse(sub$converged[i], "yes", "NO")))
      }
    }
  }
  attr(out, "B_store") <- B_store
  attr(out, "traj_store") <- traj_store
  out
}

# Plot functions

traj_to_df <- function(traj_store, field) {
  rows <- lapply(traj_store, function(r) {
    v <- r[[field]]
    if (is.null(v)) return(NULL)
    data.frame(J = r$J, Q = r$Q, seed = r$seed, method = r$method,
               floor = r$floor, iter = seq_along(v), value = v)
  })
  do.call(rbind, rows)
}

plot_log_evidence <- function(traj_store) {
  df <- traj_to_df(traj_store, "log_evidence")
  df$cell <- sprintf("J=%d Q=%d seed=%d", df$J, df$Q, df$seed)
  ggplot(df, aes(iter, value, color = method)) +
    geom_line() +
    facet_wrap(~cell, scales = "free") +
    labs(x = "EM iteration", y = "log evidence", title = "Log-evidence monotonicity")
}

plot_B_trajectory <- function(traj_store) {
  df <- traj_to_df(traj_store, "B_dist_trace")
  df$cell <- sprintf("J=%d Q=%d seed=%d", df$J, df$Q, df$seed)
  floors <- unique(df[, c("cell", "floor")])
  ggplot(df, aes(iter, value, color = method)) +
    geom_line() +
    geom_hline(data = floors, aes(yintercept = floor), linetype = "dashed", color = "grey40") +
    facet_wrap(~cell, scales = "free") +
    labs(x = "EM iteration", y = "subspace distance to B_true",
         title = "B recovery over iterations (dashed = BBP floor)")
}

plot_sigma2_trajectory <- function(traj_store, sigma2_true) {
  df <- traj_to_df(traj_store, "sigma2_trace")
  df$cell <- sprintf("J=%d Q=%d seed=%d", df$J, df$Q, df$seed)
  ggplot(df, aes(iter, value, color = method)) +
    geom_line() +
    geom_hline(yintercept = sigma2_true, linetype = "dashed", color = "grey40") +
    facet_wrap(~cell, scales = "free") +
    labs(x = "EM iteration", y = "sigma2 estimate",
         title = sprintf("sigma2 over iterations (dashed = true value %.2f)", sigma2_true))
}

plot_u_frac_trajectory <- function(traj_store) {
  df <- traj_to_df(traj_store, "u_frac_trace")
  if (is.null(df) || nrow(df) == 0) {
    message("u_frac_trace not found in traj_store; add it to fit_pfa_wbonly/fit_pfa_woodbury_lam_corr first")
    return(invisible(NULL))
  }
  df$cell <- sprintf("J=%d Q=%d seed=%d", df$J, df$Q, df$seed)
  ggplot(df, aes(iter, value, color = method)) +
    geom_line() +
    facet_wrap(~cell, scales = "free") +
    labs(x = "EM iteration", y = "u_frac",
         title = "Energy fraction along the unidentified 1_Q direction, by iteration")
}

plot_runtime <- function(res) {
  df <- res[res$method != "oracle", ]
  if (length(unique(df$J)) == 1 && length(unique(df$Q)) == 1) {
    ggplot(df, aes(method, time_s, fill = method)) +
      geom_boxplot() +
      labs(x = NULL, y = "time (s)", title = "Runtime by method")
  } else {
    df$cell <- sprintf("J=%d Q=%d", df$J, df$Q)
    ggplot(df, aes(method, time_s, fill = method)) +
      geom_col(position = "dodge") +
      facet_wrap(~cell, scales = "free_y") +
      labs(x = NULL, y = "time (s)", title = "Runtime by method and problem size")
  }
}

