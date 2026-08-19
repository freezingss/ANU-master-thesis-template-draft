source("pfa_base.R")
source("pfa_woodbury.R")
source("pfa_woodbury_traced.R")
source("basic_functions.R")
source("sim_data.R")

Q_grid <- c(30, 50)
seeds <- 1
J_fixed <- 50
K <- 2
N_per_grp <- 12
max_iter <- 100 
ll_tol <- 1e-6
s2_tol <- 1e-6

converged_ok <- function(fit) {
  le <- fit$log_evidence
  s2 <- fit$sigma2_trace         
  n <- length(le)
  if (n < 3) return(FALSE)
  ll_ok <- abs((le[n] - le[n-1]) / le[n-1]) < ll_tol
  if (is.null(s2)) return(ll_ok)
  m <- length(s2)
  s2_ok <- abs((s2[m] - s2[m-1]) / s2[m-1]) < s2_tol
  ll_ok && s2_ok
}

run_one_cell <- function(Q, seed) {
  
  dat <- simulate_pfa_data(Q = Q, K = K, J = J_fixed,
                           N_per_group = N_per_grp, seed = seed)
  
  out <- tryCatch({
    
    t_false <- system.time(
      f_false <- fit_pfa_woodbury(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                                  max_iter = max_iter, verbose = FALSE,
                                  exact_Shat = FALSE)
    )["elapsed"]
    
    t_true <- system.time(
      f_true <- fit_pfa_woodbury(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                                 max_iter = max_iter, verbose = FALSE,
                                 exact_Shat = TRUE)
    )["elapsed"]
    
    sigma2_gen <- dat$true$sigma2
    
    data.frame(
      Q = Q,
      seed = seed,
      status = "ok",
      conv_FALSE = converged_ok(f_false),
      conv_TRUE = converged_ok(f_true),
      sigma2_FALSE = f_false$sigma2,
      sigma2_TRUE = f_true$sigma2,
      sigma2_gen = sigma2_gen,
      sigma2_rel_gap = (f_false$sigma2 - f_true$sigma2) / f_true$sigma2,
      sigma2_err_FALSE = (f_false$sigma2 - sigma2_gen) / sigma2_gen,
      sigma2_err_TRUE = (f_true$sigma2  - sigma2_gen) / sigma2_gen,
      Bdist_FALSE = subspace_dist(f_false$B, dat$true$B),
      Bdist_TRUE = subspace_dist(f_true$B,  dat$true$B),
      Bdist_gap = abs(subspace_dist(f_false$B, dat$true$B) -
                               subspace_dist(f_true$B,  dat$true$B)),
      logev_gap = abs(tail(f_false$log_evidence, 1) -
                               tail(f_true$log_evidence, 1)),
      lambda_cor_gap = abs(lambda_mean_cor(f_false$lambda_hat, dat$true$lambda) -
                               lambda_mean_cor(f_true$lambda_hat,  dat$true$lambda)),
      time_FALSE = t_false,
      time_TRUE = t_true,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(Q = Q, seed = seed, status = paste("ERROR:", conditionMessage(e)),
               stringsAsFactors = FALSE)
  })
  
  out
}

grid <- expand.grid(Q = Q_grid, seed = seeds)
results_list <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  Qi <- grid$Q[i]; si <- grid$seed[i]
  cat(sprintf("[%3d/%3d] Q=%4d seed=%2d\n", i, nrow(grid), Qi, si))
  results_list[[i]] <- run_one_cell(Qi, si)
}

results_raw <- do.call(rbind, results_list)
saveRDS(results_raw, "exact_shat_fullfit_sweep_raw.rds")

n_fail <- sum(results_raw$status != "ok")
if (n_fail > 0) {
  cat(sprintf("WARNING: %d / %d cells failed, excluded from summary.\n",
              n_fail, nrow(results_raw)))
  print(results_raw[results_raw$status != "ok", c("Q", "seed", "status")])
}
results_ok <- results_raw[results_raw$status == "ok", ]

summary_list <- lapply(split(results_ok, results_ok$Q), function(df) {
  data.frame(
    Q = df$Q[1],
    n_seeds = nrow(df),
    n_not_converged = sum(!df$conv_FALSE | !df$conv_TRUE),
    mean_sigma2_rel_gap = mean(df$sigma2_rel_gap),
    sd_sigma2_rel_gap = sd(df$sigma2_rel_gap),
    mean_sigma2_err_FALSE = mean(df$sigma2_err_FALSE),
    mean_sigma2_err_TRUE = mean(df$sigma2_err_TRUE),
    mean_Bdist_gap = mean(df$Bdist_gap),
    sd_Bdist_gap = sd(df$Bdist_gap),
    mean_logev_gap = mean(df$logev_gap),
    mean_lambda_cor_gap = mean(df$lambda_cor_gap),
    mean_time_FALSE = mean(df$time_FALSE),
    mean_time_TRUE = mean(df$time_TRUE),
    speedup = mean(df$time_TRUE) / mean(df$time_FALSE)
  )
})

results_summary <- do.call(rbind, summary_list)
results_summary <- results_summary[order(results_summary$Q), ]
rownames(results_summary) <- NULL
saveRDS(results_summary, "exact_shat_fullfit_sweep_summary.rds")
print(results_summary)

write_latex_table <- function(df, file) {
  con <- file(file, "w")
  on.exit(close(con))
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    writeLines(sprintf(
      "%d & %.4f (%.4f) & %.4f (%.4f) & %.4f & %.4f & %.1f$\\times$ \\\\",
      r$Q,
      r$mean_sigma2_rel_gap, r$sd_sigma2_rel_gap,
      r$mean_Bdist_gap,      r$sd_Bdist_gap,
      r$mean_logev_gap,
      r$mean_lambda_cor_gap,
      r$speedup
    ), con)
  }
}

dir.create("tables", showWarnings = FALSE)
write_latex_table(results_summary, "tables/exact_shat_fullfit_gap_by_Q.tex")