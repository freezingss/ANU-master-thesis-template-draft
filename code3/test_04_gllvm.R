library(gllvm)
library(bench) # timing

source("pfa_base_multinomial.R") # base multinomial
source("pfa_woodbury_only.R") # woodbury
source("pfa_woodbury_lam_corrected.R") # woodbury + lambda correction
source("sim_data.R")
source("basic_functions.R")

if (!exists("fit_gllvm_ref")) {
  
  fit_gllvm_ref <- function(Y, X, group, K, verbose = FALSE) {
    if (!requireNamespace("gllvm", quietly = TRUE)) {
      return(list(B = NULL, sigma2 = NA_real_, time = NA_real_,
                  ok = FALSE, converged = FALSE, error = "gllvm package not installed"))
    }
    t0 <- proc.time()["elapsed"]
    
    ug <- sort(unique(group))
    Yg <- t(sapply(ug, function(g) colSums(Y[group == g, , drop = FALSE])))
    Xg <- NULL
    if (!is.null(X)) {
      Xm <- as.matrix(X)
      Xg <- t(sapply(ug, function(g) colMeans(Xm[group == g, , drop = FALSE])))
    }
    
    fit <- tryCatch(
      gllvm::gllvm(y = Yg, X = Xg, family = poisson(), num.lv = K,
                   method = "LA", starting.val = "res"),
      error = function(e) e)
    
    el <- proc.time()["elapsed"] - t0
    if (inherits(fit, "error")) {
      return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE,
                  converged = FALSE, error = conditionMessage(fit)))
    }
    
    B <- tryCatch(as.matrix(gllvm::getLoadings(fit)), error = function(e) NULL)
    if (is.null(B)) {
      B <- tryCatch(as.matrix(fit$params$theta), error = function(e) NULL)
    }
    if (!is.null(B) && nrow(B) != ncol(Yg)) B <- t(B) 
    
    conv <- tryCatch(is.finite(fit$logL) && (is.null(fit$convergence) || isTRUE(fit$convergence == 0)),
                     error = function(e) FALSE)
    
    if (is.null(B)) {
      return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE,
                  converged = conv, error = "loading extraction failed"))
    }
    list(B = B, sigma2 = NA_real_, time = el, ok = TRUE, converged = conv,
         error = NA_character_)
  }
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
Js <- c(25, 50, 100, 150)
Qs <- c(50, 100, 150, 200, 250)
seeds <- 1:3
K <- 2
Nj <- 15

gllvm_results <- data.frame(
  J = integer(), Q = integer(), seed = integer(), method = character(),
  time_s = numeric(), swigma2 = numeric(), d_true = numeric(),
  converged = logical(), log_lik_monotone = logical(),
  stringsAsFactors = FALSE
)

cat(sprintf("%5s %6s %6s | %14s %10s %10s %10s %8s\n",
            "J", "Q", "seed", "method", "time(s)", "sigma2", "d vs true", "conv?"))

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
      # 
      # b_wb <- bench::mark(
      #   fw_only <<- fit_pfa_wbonly_traced(dat$Y, dat$X, dat$group, K = K, M = dat$M,
      #                                     max_iter = 80, tol = 1e-3, sigma2_init = target, 
      #                                     verbose = FALSE, estep_max_iter = 100,
      #                                     estep_gtol = 1e-3, B_true = NULL),
      #   iterations = 1, check = FALSE, memory = FALSE)
      # t_wb <- as.numeric(b_wb$median)
      
      # REML-corrected sigma2
      # s2_reml <- fw_old$sigma2 * J / (J - K)
      
      # corrected (Woodbury + Edgeworth)
      # b_corr <- bench::mark(
      #   fw_corr <<- fit_pfa_woodbury_lam_corr(dat$Y, dat$X, dat$group, K, M = dat$M,
      #                                         max_iter = 80, tol = 1e-3, sigma2_init = target,
      #                                         verbose = FALSE, estep_max_iter = 100, estep_gtol = 1e-3,
      #                                         B_true = NULL,
      #                                         corr_max_rel = 0.5, 
      #                                         use_lambda_correction = TRUE),
      #   iterations = 1, check = FALSE, memory = FALSE)
      # t_corr <- as.numeric(b_corr$median)
      
      # gllvm()
      b_gl <- tryCatch(
        bench::mark(fr <<- fit_gllvm_ref(dat$Y, dat$X, dat$group, K = K, verbose = FALSE),
                    iterations = 1, check = FALSE, memory = FALSE),
        error = function(e) NULL)
      ok <- !is.null(b_gl) && isTRUE(fr$ok)
      t_gl <- if (ok) as.numeric(b_gl$median) else NA
      
      # d_base_true <- subspace_dist(fw_base$B, dat$true$B)
      # d_only_true  <- subspace_dist(fw_only$B, dat$true$B)
      # d_corr_true <- subspace_dist(fw_corr$B, dat$true$B)
      # d_G_true    <- subspace_dist(fw_G$B, dat$true$B)
      # d_Gf_true   <- subspace_dist(fw_G_frz$B, dat$true$B)
      d_gl_true   <- if (ok) subspace_dist(fr$B, dat$true$B) else NA
      
      # rows <- data.frame(
      #   J = J, Q = Q, seed = seed,
      #   method = c("base", "woodbury only", "corrected", "gllvm()"),
      #   time_s = c(NA, t_wb, t_corr, t_gl),
      #   sigma2 = c(NA, fw_only$sigma2, fw_corr$sigma2, if (ok) fr$sigma2 else NA),
      #   d_true = c(NA, d_only_true, d_corr_true, d_gl_true),
      #   converged = c(NA, fw_only$converged, fw_corr$converged,
      #                 if (ok) fr$converged else NA),
      #   log_lik_monotone = c(sapply(list(NA, fw_old, fw_corr), function(f) monotone_diagnostics(f)$monotone), NA),
      #   stringsAsFactors = FALSE
      # )
      
      rows <- data.frame(
        J = J, Q = Q, seed = seed,
        method = c("gllvm()"),
        time_s = c(t_gl),
        sigma2 = c(if (ok) fr$sigma2 else NA),
        d_true = c(d_gl_true),
        converged = c(if (ok) fr$converged else NA),
        stringsAsFactors = FALSE
      )
      q_results <- rbind(gllvm_results, rows)
      
      for (i in seq_len(nrow(rows))) {
        cat(sprintf("%5d %6d %6d | %14s %10.2f %10.4f %10s %8s\n",
                    J, Q, seed, rows$method[i], rows$time_s[i], rows$sigma2[i],
                    if (is.na(rows$d_true[i])) "-" else sprintf("%.4f", rows$d_true[i]),
                    if (is.na(rows$converged[i])) "-" else ifelse(rows$converged[i], "yes", "NO")))
      }
      cat(sprintf("  J=%d Q=%d seed=%d | true sigma2=%.4f\n", J, Q, seed, dat$true$sigma2))
    }
  }
}

# # Calculate the mean of estimations
# avg_time  <- aggregate(time_s ~ method + Q + J, data = q_results,
#                        FUN = function(x) mean(x, na.rm = TRUE))
# avg_sigma <- aggregate(sigma2 ~ method + Q + J, data = q_results,
#                        FUN = function(x) mean(x, na.rm = TRUE))
# avg_d     <- aggregate(d_true ~ method + Q + J, data = q_results,
#                        FUN = function(x) mean(x, na.rm = TRUE))
# 
# avg_results <- Reduce(function(a, b) merge(a, b, by = c("method", "Q", "J")),
#                       list(avg_time, avg_sigma, avg_d))
# avg_results <- avg_results[order(avg_results$Q, avg_results$J, avg_results$method), ]
# 
# cat("\n Averaged across seeds, per (method, Q, J) \n")
# cat(sprintf("%14s %6s %6s | %10s %10s %10s\n",
#             "method", "Q", "J", "time(s)", "sigma2", "d_true"))
# for (i in seq_len(nrow(avg_results))) {
#   cat(sprintf("%14s %6d %6d | %10.2f %10.4f %10.4f\n",
#               avg_results$method[i], avg_results$Q[i], avg_results$J[i],
#               avg_results$time_s[i], avg_results$sigma2[i], avg_results$d_true[i]))
# }

