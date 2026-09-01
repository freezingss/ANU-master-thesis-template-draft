library(glmmTMB)
library(bench)

stopifnot(packageVersion("glmmTMB") >= "1.1.8")

source("pfa_base_multinomial.R") # base model
source("pfa_woodbury_only.R") # woodbury acceleration
source("pfa_woodbury_lam_corrected.R") # woodbury + lambda correction
source("sim_data.R")
source("basic_functions.R")

fit_glmmTMB_ref <- function(Y, X, group, K, time_limit = Inf, verbose = FALSE) {
  long <- build_long(Y, X, group)
  t0 <- proc.time()["elapsed"]
  
  n_obs_levels <- length(unique(long$obs))
  n_gc_pairs <- length(unique(interaction(long$group, long$category, drop = TRUE)))
  if (n_obs_levels != n_gc_pairs) {
    stop(sprintf("obs is not unique per (group, category): %d obs levels vs %d (group, category) pairs.",
                 n_obs_levels, n_gc_pairs))
  }
  if (anyNA(long$count) || anyNA(long$category) || anyNA(long$group) ||
      anyNA(long$obs) || anyNA(long$log_total)) stop("NA found in build_long() output.")
  if (any(!is.finite(long$log_total))) stop("long$log_total has non-finite value(s).")
  
  form <- as.formula(paste0("count ~ category + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))
  fit_once <- function(ctrl) {
    tryCatch(glmmTMB(form, data = long, family = poisson(link = "log"),
                     offset = long$log_total, control = ctrl),
             error = function(e) e)
  }
  ctrl_default <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  fit <- fit_once(ctrl_default)
  if (inherits(fit, "error")) {
    ctrl_res <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000),
                               start_method = list(method = "res"))
    fit2 <- fit_once(ctrl_res)
    if (inherits(fit2, "error")) {
      el <- proc.time()["elapsed"] - t0
      return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE, converged = FALSE,
                  error = paste0("default: ", conditionMessage(fit), " | res-start: ", conditionMessage(fit2))))
    }
    fit <- fit2
  }
  el <- proc.time()["elapsed"] - t0
  conv <- isTRUE(fit$sdr$pdHess)
  
  L <- tryCatch({
    vc <- as.matrix(VarCorr(fit)$cond$group)
    e  <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
    e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(e$values[1:K], 0)), K)
  }, error = function(e) NULL)
  if (is.null(L)) {
    L <- tryCatch(as.matrix(fit$obj$env$report(fit$fit$parfull)$fact_load[[1]]),
                  error = function(e) NULL)
  }
  if (is.null(L)) {
    return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE, converged = conv,
                error = "loading extraction failed"))
  }
  sigma2_hat <- tryCatch(as.numeric(VarCorr(fit)$cond$obs)[1], error = function(e) NA_real_)
  list(B = L, sigma2 = sigma2_hat, time = el, ok = TRUE, converged = conv, error = NA_character_)
}

fit_pfa_woodbury_G_wrapper <- function(Y, X, group, K, M, max_iter = 80, fix_sigma2 = NULL) {
  fit_pfa_woodbury_G(Y, X, group, K, M = M, max_iter = max_iter, tol = 1e-8,
                      sigma2_init = 0.3, verbose = FALSE, fix_sigma2 = fix_sigma2)
}

# Monotone Check
is_monotone_loglik <- function(fit) {
  all(diff(fit$log_evidence) >= 0)
}

monotone_diagnostics <- function(fit) {
  d <- diff(fit$log_evidence)
  list(monotone = all(d >= 0), n_decreasing = sum(d < 0),
       first_decrease_iter = if (any(d < 0)) which(d < 0)[1] + 1 else NA_integer_)
}

# RUN CODE
Js    <- c(100, 150)
Qs    <- c(50, 100, 150)
# Qs <- c(50)
seeds <- 1:3
# seeds <- c(1,3) # Observed severe unstable for old woodbury + PME
K <- 2
Nj <- 15

q_results <- data.frame(
  J = integer(), Q = integer(), seed = integer(), method = character(),
  time_s = numeric(), sigma2 = numeric(), d_true = numeric(), converged = logical(),
  log_lik_monotone = logical(),
  stringsAsFactors = FALSE
)

cat(sprintf("%5s %6s %6s | %14s %10s %10s %10s %8s\n %10s\n",
            "J", "Q", "seed", "method", "time(s)", "sigma2", "d vs true", "conv?", "LL mono?"))

for (J in Js) {
  for (Q in Qs) {
    for (seed in seeds) {
      
      # target = (J - K) / J # Large sigma2 test
      target = 0.3 # Small sigma2 test

      dat <- simulate_pfa_data(Q = Q, K = K, J = J, P = 3, N_per_group = Nj, sigma2 = target, M_rate = 150, seed = seed)
      
      # b_base <- bench::mark(
      #   fw_base <<- fit_pfa_dense(dat$Y, dat$X, dat$group, K = K, M = dat$M,
      #                             max_iter = 80, tol = 1e-3, sigma2_init = target, 
      #                             verbose = FALSE, estep_max_iter = 100,
      #                             estep_gtol = 1e-3),
      #   iterations = 1, check = FALSE, memory = FALSE)
      # t_base <-as.numeric(b_base$median)

      b_wb <- bench::mark(
        fw_only <<- fit_pfa_wbonly_traced(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                                     max_iter = 80, tol = 1e-3, sigma2_init = target, 
                                     verbose = FALSE, estep_max_iter = 100,
                                     estep_gtol = 1e-3, B_true = NULL),
        iterations = 1, check = FALSE, memory = FALSE)
      t_wb <- as.numeric(b_wb$median)

      # REML-corrected sigma2
      # s2_reml <- fw_old$sigma2 * J / (J - K)

      # corrected (Woodbury + Edgeworth)
      b_corr <- bench::mark(
        fw_corr <<- fit_pfa_woodbury_lam_corr(dat$Y, dat$X, dat$group, K, M = dat$M,
                                               max_iter = 80, tol = 1e-3, sigma2_init = target,
                                               verbose = FALSE, estep_max_iter = 100, estep_gtol = 1e-3,
                                               B_true = NULL,
                                               corr_max_rel = 0.5, 
                                               use_lambda_correction = TRUE),
        iterations = 1, check = FALSE, memory = FALSE)
      t_corr <- as.numeric(b_corr$median)

      # # new (G, unfrozen)
      # b_G <- bench::mark(
      #   fw_G <<- fit_pfa_woodbury_G_wrapper(dat$Y, dat$X, dat$group, K, dat$M,
      #                                        max_iter = 80, fix_sigma2 = NULL),
      #   iterations = 1, check = FALSE, memory = FALSE)
      # t_G <- as.numeric(b_G$median)
      # 
      # # new (G, frozen at REML value)
      # b_G_frz <- bench::mark(
      #   fw_G_frz <<- fit_pfa_woodbury_G_wrapper(dat$Y, dat$X, dat$group, K, dat$M,
      #                                            max_iter = 80, fix_sigma2 = s2_reml),
      #   iterations = 1, check = FALSE, memory = FALSE)
      # t_G_frz <- as.numeric(b_G_frz$median)

      # glmmTMB rr()
      b_gt <- tryCatch(
        bench::mark(fr <<- fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K = K, verbose = FALSE),
                    iterations = 1, check = FALSE, memory = FALSE),
        error = function(e) NULL)
      ok <- !is.null(b_gt) && isTRUE(fr$ok)
      t_gt <- if (ok) as.numeric(b_gt$median) else NA
      
      # d_base_true <- subspace_dist(fw_base$B, dat$true$B)
      d_only_true  <- subspace_dist(fw_only$B, dat$true$B)
      d_corr_true <- subspace_dist(fw_corr$B, dat$true$B)
      # d_G_true    <- subspace_dist(fw_G$B, dat$true$B)
      # d_Gf_true   <- subspace_dist(fw_G_frz$B, dat$true$B)
      d_gt_true   <- if (ok) subspace_dist(fr$B, dat$true$B) else NA

      # rows <- data.frame(
      #   J = J, Q = Q, seed = seed,
      #   method = c("old (Woodbury EM)", "corrected (Woodbury+Edgeworth)",
      #              "new (G, unfrozen)", "new (G, frozen)", "glmmTMB rr()"),
      #   time_s = c(t_wb, t_corr, t_G, t_G_frz, t_gt),
      #   d_true = c(d_old_true, d_corr_true, d_G_true, d_Gf_true, d_gt_true),
      #   d_vs_wb_old = c(0,
      #                    subspace_dist(fw_corr$B, fw_old$B),
      #                    subspace_dist(fw_G$B, fw_old$B),
      #                    subspace_dist(fw_G_frz$B, fw_old$B),
      #                    if (ok) subspace_dist(fr$B, fw_old$B) else NA),
      #   converged = c("-", "-", "-", "-", if (ok) ifelse(isTRUE(fr$converged), "yes", "NO") else "FAILED"),
      #   stringsAsFactors = FALSE
      # )
      rows <- data.frame(
        J = J, Q = Q, seed = seed,
        method = c("woodbury only", "corrected", "glmmTMB rr()"),
        time_s = c(t_wb, t_corr, t_gt),
        sigma2 = c(fw_only$sigma2, fw_corr$sigma2, if (ok) fr$sigma2 else NA),
        d_true = c(d_only_true, d_corr_true, d_gt_true),
        converged = c(fw_only$converged, fw_corr$converged,
                      if (ok) fr$converged else NA),
        log_lik_monotone = c(sapply(list(fw_only, fw_corr), function(f) monotone_diagnostics(f)$monotone), NA),
        stringsAsFactors = FALSE
      )
      q_results <- rbind(q_results, rows)
      
      for (i in seq_len(nrow(rows))) {
        cat(sprintf("%5d %6d %6d | %14s %10.2f %10.4f %10s %8s\n %10s\n",
                    J, Q, seed, rows$method[i], rows$time_s[i], rows$sigma2[i],
                    if (is.na(rows$d_true[i])) "-" else sprintf("%.4f", rows$d_true[i]),
                    if (is.na(rows$converged[i])) "-" else ifelse(rows$converged[i], "yes", "NO"),
                    if (is.na(rows$log_lik_monotone[i])) "-" else ifelse(rows$log_lik_monotone[i], "yes", "NO")))
      }
      cat(sprintf("  J=%d Q=%d seed=%d | true sigma2=%.4f\n", J, Q, seed, dat$true$sigma2))
    }
  }
}

# Calculate the mean of estimations
avg_time  <- aggregate(time_s ~ method + Q + J, data = q_results,
                       FUN = function(x) mean(x, na.rm = TRUE))
avg_sigma <- aggregate(sigma2 ~ method + Q + J, data = q_results,
                       FUN = function(x) mean(x, na.rm = TRUE))
avg_d     <- aggregate(d_true ~ method + Q + J, data = q_results,
                       FUN = function(x) mean(x, na.rm = TRUE))

avg_results <- Reduce(function(a, b) merge(a, b, by = c("method", "Q", "J")),
                      list(avg_time, avg_sigma, avg_d))
avg_results <- avg_results[order(avg_results$Q, avg_results$J, avg_results$method), ]

cat("\n Averaged across seeds, per (method, Q, J) \n")
cat(sprintf("%14s %6s %6s | %10s %10s %10s\n",
            "method", "Q", "J", "time(s)", "sigma2", "d_true"))
for (i in seq_len(nrow(avg_results))) {
  cat(sprintf("%14s %6d %6d | %10.2f %10.4f %10.4f\n",
              avg_results$method[i], avg_results$Q[i], avg_results$J[i],
              avg_results$time_s[i], avg_results$sigma2[i], avg_results$d_true[i]))
}

