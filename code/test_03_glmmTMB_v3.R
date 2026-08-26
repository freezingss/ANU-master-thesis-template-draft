library(glmmTMB)
library(bench)

stopifnot(packageVersion("glmmTMB") >= "1.1.8")

source("pfa_woodbury.R") # fit_pfa_woodbury()
# source("pfa_woodbury_G.R") # fit_pfa_woodbury_G()
# source("pfa_woodbury_lam_corr.R") # fit_pfa_woodbury_lam_corr()
# source("pfa_woodbury_lambda_corrected.R") # correct_lambda_edgeworth()
source("pfa_woodbury_G_dense_final.R")
source("pfa_woodbury_lam_corrected_final.R")
source("sim_data.R")
source("basic_functions.R")

# if (!exists("fit_glmmTMB_ref")) {
#   fit_glmmTMB_ref <- function(Y, X, group, K, time_limit = Inf, verbose = TRUE) {
#     long <- build_long(Y, X, group)
#     t0 <- proc.time()["elapsed"]
# 
#     n_obs_levels <- length(unique(long$obs))
#     n_gc_pairs   <- length(unique(interaction(long$group, long$category, drop = TRUE)))
#     if (n_obs_levels != n_gc_pairs) {
#       stop(sprintf("obs is not unique per (group, category): %d obs levels vs %d (group, category) pairs.",
#                     n_obs_levels, n_gc_pairs))
#     }
#     if (anyNA(long$count) || anyNA(long$category) || anyNA(long$group) ||
#         anyNA(long$obs) || anyNA(long$log_total)) stop("NA found in build_long() output.")
#     if (any(!is.finite(long$log_total))) stop("long$log_total has non-finite value(s).")
# 
#     form <- as.formula(paste0("count ~ category + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))
#     fit_once <- function(ctrl) {
#       tryCatch(glmmTMB(form, data = long, family = poisson(link = "log"),
#                         offset = long$log_total, control = ctrl),
#                 error = function(e) e)
#     }
#     ctrl_default <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
#     fit <- fit_once(ctrl_default)
#     if (inherits(fit, "error")) {
#       ctrl_res <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000),
#                                   start_method = list(method = "res"))
#       fit2 <- fit_once(ctrl_res)
#       if (inherits(fit2, "error")) {
#         el <- proc.time()["elapsed"] - t0
#         return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE, converged = FALSE,
#                     error = paste0("default: ", conditionMessage(fit), " | res-start: ", conditionMessage(fit2))))
#       }
#       fit <- fit2
#     }
#     el <- proc.time()["elapsed"] - t0
#     conv <- isTRUE(fit$sdr$pdHess)
# 
#     L <- tryCatch({
#       vc <- as.matrix(VarCorr(fit)$cond$group)
#       e  <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
#       e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(e$values[1:K], 0)), K)
#     }, error = function(e) NULL)
#     if (is.null(L)) {
#       L <- tryCatch(as.matrix(fit$obj$env$report(fit$fit$parfull)$fact_load[[1]]),
#                      error = function(e) NULL)
#     }
#     if (is.null(L)) {
#       return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE, converged = conv,
#                   error = "loading extraction failed"))
#     }
#     sigma2_hat <- tryCatch(as.numeric(VarCorr(fit)$cond$obs)[1], error = function(e) NA_real_)
#     list(B = L, sigma2 = sigma2_hat, time = el, ok = TRUE, converged = conv, error = NA_character_)
#   }
# }

fit_pfa_woodbury_G_wrapper <- function(Y, X, group, K, M, max_iter = 80, fix_sigma2 = NULL) {
  fit_pfa_woodbury_G(Y, X, group, K, M = M, max_iter = max_iter, tol = 1e-8,
                      sigma2_init = 0.3, verbose = FALSE, fix_sigma2 = fix_sigma2)
}

# RUN CODE
Js    <- c(50, 100)
# Qs    <- c(25, 50, 75, 100)
Qs <- c(50)
# seeds <- 1:3
seeds <- c(1,3)
K <- 2
Nj <- 15

q_results <- data.frame(
  J = integer(), Q = integer(), seed = integer(), method = character(),
  time_s = numeric(),
  d_true = numeric(), d_vs_wb_old = numeric(), converged = character(),
  stringsAsFactors = FALSE
)

cat(sprintf("%5s %6s %6s | %30s %14s %10s %14s %8s\n",
            "J", "Q", "seed", "method", "time(s)", "d vs true", "d vs wb-old", "conv?"))

for (J in Js) {
  for (Q in Qs) {
    for (seed in seeds) {
      
      target = (J - K) / J

      dat <- simulate_pfa_data(Q = Q, K = K, J = J, N_per_group = Nj, sigma2 = target, seed = seed)

      # old (Woodbury EM) unmodified baseline
      b_wb <- bench::mark(
        fw_old <<- fit_pfa_woodbury(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                                     max_iter = 80, verbose = FALSE, exact_Shat = FALSE),
        iterations = 1, check = FALSE, memory = FALSE)
      t_wb <- as.numeric(b_wb$median)

      # REML-corrected sigma2
      s2_reml <- fw_old$sigma2 * J / (J - K)

      # corrected (Woodbury + Edgeworth) unfrozen sigma2
      b_corr <- bench::mark(
        fw_corr <<- fit_pfa_woodbury_lam_corr(dat$Y, dat$X, dat$group, K, M = dat$M,
                                               max_iter = 80, tol = 1e-8, sigma2_init = 0.3,
                                               verbose = FALSE, exact_Shat = FALSE,
                                               use_lambda_correction = TRUE),
        iterations = 1, check = FALSE, memory = FALSE)
      t_corr <- as.numeric(b_corr$median)

      # new (G, unfrozen)
      b_G <- bench::mark(
        fw_G <<- fit_pfa_woodbury_G_wrapper(dat$Y, dat$X, dat$group, K, dat$M,
                                             max_iter = 80, fix_sigma2 = NULL),
        iterations = 1, check = FALSE, memory = FALSE)
      t_G <- as.numeric(b_G$median)

      # new (G, frozen at REML value)
      b_G_frz <- bench::mark(
        fw_G_frz <<- fit_pfa_woodbury_G_wrapper(dat$Y, dat$X, dat$group, K, dat$M,
                                                 max_iter = 80, fix_sigma2 = s2_reml),
        iterations = 1, check = FALSE, memory = FALSE)
      t_G_frz <- as.numeric(b_G_frz$median)

      # glmmTMB rr()
      b_gt <- tryCatch(
        bench::mark(fr <<- fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K = K, verbose = FALSE),
                    iterations = 1, check = FALSE, memory = FALSE),
        error = function(e) NULL)
      ok <- !is.null(b_gt) && isTRUE(fr$ok)
      t_gt <- if (ok) as.numeric(b_gt$median) else NA

      d_old_true  <- subspace_dist(fw_old$B, dat$true$B)
      d_corr_true <- subspace_dist(fw_corr$B, dat$true$B)
      d_G_true    <- subspace_dist(fw_G$B, dat$true$B)
      d_Gf_true   <- subspace_dist(fw_G_frz$B, dat$true$B)
      d_gt_true   <- if (ok) subspace_dist(fr$B, dat$true$B) else NA

      rows <- data.frame(
        J = J, Q = Q, seed = seed,
        method = c("old (Woodbury EM)", "corrected (Woodbury+Edgeworth)",
                   "new (G, unfrozen)", "new (G, frozen)", "glmmTMB rr()"),
        time_s = c(t_wb, t_corr, t_G, t_G_frz, t_gt),
        d_true = c(d_old_true, d_corr_true, d_G_true, d_Gf_true, d_gt_true),
        d_vs_wb_old = c(0,
                         subspace_dist(fw_corr$B, fw_old$B),
                         subspace_dist(fw_G$B, fw_old$B),
                         subspace_dist(fw_G_frz$B, fw_old$B),
                         if (ok) subspace_dist(fr$B, fw_old$B) else NA),
        converged = c("-", "-", "-", "-", if (ok) ifelse(isTRUE(fr$converged), "yes", "NO") else "FAILED"),
        stringsAsFactors = FALSE
      )
      q_results <- rbind(q_results, rows)

      for (i in seq_len(nrow(rows))) {
        cat(sprintf("%5d %6d %6d | %30s %14.2f %10s %14s %8s\n",
                    J, Q, seed, rows$method[i], rows$time_s[i],
                    if (is.na(rows$d_true[i])) "-" else sprintf("%.4f", rows$d_true[i]),
                    if (is.na(rows$d_vs_wb_old[i])) "-" else sprintf("%.4f", rows$d_vs_wb_old[i]),
                    rows$converged[i]))
      }
      cat(sprintf("  J=%d Q=%d seed=%d | wb s2=%.4f  corr s2=%.4f  G s2=%.4f  G+frz s2=%.4f (s2_reml=%.4f)  true s2=%.4f\n",
                  J, Q, seed, fw_old$sigma2, fw_corr$sigma2, fw_G$sigma2, fw_G_frz$sigma2, s2_reml, dat$true$sigma2))
    }
  }
}

saveRDS(q_results, "q_results_J50100(50)_Q50_seed13_largesigma2.rds")

source("diagnose_sigma.R")

# # Diagnostic analysis
# 
# source("diagnose_group.R")
# 
# inspect_group_data(20, dat$Y, dat$group, dat$M)
# 
# tier5_rate <- apply(fw_G$tier_used == 5, 2, mean)   
# hhi_check  <- inspect_group_data(1, dat$Y, dat$group, dat$M)$all_hhi
# 
# cor(tier5_rate, hhi_check)
# plot(hhi_check, tier5_rate, xlab = "HHI (category concentration)", ylab = "tier-5 fallback rate",
#      main = "concentration vs fallback rate, all groups")
# text(hhi_check[20], tier5_rate[20], "20", col = "red", pos = 3)
# 
# diagnose_group_trace(20, dat$Y, dat$X, dat$group, K, M = dat$M,
#                      max_iter = 80, tol = 1e-10, fix_sigma2 = NULL,
#                      B_true = dat$true$B)
# 
# source("diagnose_theta_breakdown.R")
# 
# diagnose_theta_breakdown(20, dat$Y, dat$X, dat$group, K, M = dat$M,
#                          max_iter = 80, tol = 1e-10, fix_sigma2 = NULL,
#                          B_true = dat$true$B,
#                          window_start = 40, window_end = 65)