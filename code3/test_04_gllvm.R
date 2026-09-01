# gllvm() result: compared with wb, corr and glmmTMB() in test_03_glmmTMB  

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
      # gllvm内部会用列名去拼公式字符串, Xg没有列名的话拼出来的公式会变成
      # "~ 0 + "这种缺变量名的残缺字符串, parse时直接报unexpected end of
      # input。这里强制给它编上列名, 不管Xm本身有没有列名。
      colnames(Xg) <- if (!is.null(colnames(Xm))) colnames(Xm) else paste0("x", seq_len(ncol(Xg)))
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
    # sigma2 这里恒为 NA: gllvm 的 Poisson/NB 潜变量模型结构上没有一个跟
    # FA-DMR 的 sigma2*I_Q (idiosyncratic 噪声) 直接对应的量, 不去凑一个假的。
    list(B = B, sigma2 = NA_real_, time = el, ok = TRUE, converged = conv,
         error = NA_character_)
  }
}

# wb/corr 才需要看log_evidence是否单调, gllvm的拟合不是EM迭代, 没有这个概念,
# 这两个函数在这个文件里已经用不到了, 一起注释掉。
# is_monotone_loglik <- function(fit) {
#   all(diff(fit$log_evidence) >= 0)
# }
#
# monotone_diagnostics <- function(fit) {
#   d <- diff(fit$log_evidence)
#   list(monotone = all(d >= 0), n_decreasing = sum(d < 0),
#        first_decrease_iter = if (any(d < 0)) which(d < 0)[1] + 1 else NA_integer_)
# }

# RUN CODE
Js <- c(25, 50, 100, 150)
Qs <- c(50, 100, 150, 200, 250)
seeds <- 1:3
K <- 2
Nj <- 15

gllvm_results <- data.frame(
  J = integer(), Q = integer(), seed = integer(),
  method = character(), time_s = numeric(),
  sigma2 = numeric(), d_true = numeric(), converged = logical(),
  error_msg = character(), stringsAsFactors = FALSE
)

cat(sprintf("%5s %6s %6s | %10s %10s %10s %8s\n",
            "J", "Q", "seed", "time(s)", "sigma2", "d vs true", "conv?"))

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
      
      # wb (woodbury only) —— 已经在glmmTMB那份脚本里跑过, 这里不重复跑
      # b_wb <- bench::mark(
      #   fw_only <<- fit_pfa_wbonly_traced(dat$Y, dat$X, dat$group, K = K, M = dat$M,
      #                                     max_iter = 80, tol = 1e-3, sigma2_init = target, 
      #                                     verbose = FALSE, estep_max_iter = 100,
      #                                     estep_gtol = 1e-3, B_true = NULL),
      #   iterations = 1, check = FALSE, memory = FALSE)
      # t_wb <- as.numeric(b_wb$median)
      
      # REML-corrected sigma2
      # s2_reml <- fw_old$sigma2 * J / (J - K)
      
      # corrected (Woodbury + Edgeworth) —— 已经在glmmTMB那份脚本里跑过, 这里不重复跑
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
      d_gl_true   <- if (ok) subspace_dist(fr$B, dat$true$B) else NA
      
      row <- data.frame(
        J = J, Q = Q, seed = seed,
        method = "gllvm()",
        time_s = t_gl,
        sigma2 = if (ok) fr$sigma2 else NA_real_,
        d_true = d_gl_true,
        converged = if (ok) fr$converged else NA,
        error_msg = if (!is.null(b_gl)) fr$error else "bench::mark itself errored",
        stringsAsFactors = FALSE
      )
      gllvm_results <- rbind(gllvm_results, row)
      
      cat(sprintf("%5d %6d %6d | %10.2f %10s %10s %8s\n",
                  J, Q, seed, row$time_s,
                  if (is.na(row$sigma2)) "-" else sprintf("%.4f", row$sigma2),
                  if (is.na(row$d_true)) "-" else sprintf("%.4f", row$d_true),
                  if (is.na(row$converged)) "-" else ifelse(row$converged, "yes", "NO")))
      if (!ok) cat(sprintf("    -> gllvm error: %s\n", row$error_msg))
      cat(sprintf("  J=%d Q=%d seed=%d | true sigma2=%.4f\n", J, Q, seed, dat$true$sigma2))
    }
  }
}

# # Calculate the mean of estimations
# avg_time  <- aggregate(time_s ~ Q + J, data = gllvm_results,
#                        FUN = function(x) mean(x, na.rm = TRUE))
# avg_d     <- aggregate(d_true ~ Q + J, data = gllvm_results,
#                        FUN = function(x) mean(x, na.rm = TRUE))
#
# avg_results <- merge(avg_time, avg_d, by = c("Q", "J"))
# avg_results <- avg_results[order(avg_results$Q, avg_results$J), ]
#
# cat("\n Averaged across seeds, per (Q, J) \n")
# cat(sprintf("%6s %6s | %10s %10s\n", "Q", "J", "time(s)", "d_true"))
# for (i in seq_len(nrow(avg_results))) {
#   cat(sprintf("%6d %6d | %10.2f %10.4f\n",
#               avg_results$Q[i], avg_results$J[i],
#               avg_results$time_s[i], avg_results$d_true[i]))
# }



# Test
dat <- simulate_pfa_data(Q = 50, K = 2, J = 25, P = 3, N_per_group = 15, sigma2 = 0.3, M_rate = 150, seed = 1)
is.null(dat$X)
dim(dat$Y)
colnames(dat$Y)

set.seed(1)
Ytoy <- matrix(rpois(200, lambda = 5), nrow = 20, ncol = 10)
gllvm::gllvm(y = Ytoy, family = poisson(), num.lv = 2, method = "LA")

gllvm::gllvm(y = dat$Y, family = poisson(), num.lv = 2, method = "LA")
traceback()

# Test 2
ug <- sort(unique(dat$group))
Yg <- t(sapply(ug, function(g) colSums(dat$Y[dat$group == g, , drop = FALSE])))
Xm <- as.matrix(dat$X)
Xg <- t(sapply(ug, function(g) colMeans(Xm[dat$group == g, , drop = FALSE])))
colnames(Xg) <- if (!is.null(colnames(Xm))) colnames(Xm) else paste0("x", seq_len(ncol(Xg)))

dim(Xg)
colnames(Xg)

gllvm::gllvm(y = Yg, X = Xg, family = poisson(), num.lv = 2, method = "LA")
traceback()
