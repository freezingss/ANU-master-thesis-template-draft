library(gllvm)
library(bench)

source("pfa_woodbury.R")
source("sim_data.R")
source("basic_functions.R")


fit_gllvm_ref <- function(Y, X, group, K) {
  Yg <- rowsum(Y, group)               
  t0 <- proc.time()["elapsed"]
  fit <- tryCatch(
    gllvm(y = Yg, family = poisson(), num.lv = K, method = "LA",
          control = list(reltol = 1e-8, maxit = 1000)),
    error = function(e) NULL)
  el <- proc.time()["elapsed"] - t0
  if (is.null(fit)) return(list(B = NULL, time = el, ok = FALSE, converged = FALSE))
  list(B = fit$params$theta, time = el, ok = TRUE,
       converged = isTRUE(fit$convergence == 0))
}

cat(sprintf("%6s | %14s %14s %14s %14s %10s %10s %10s %8s\n",
            "Q", "Woodbury(s)", "wb mem(MB)", "gllvm(s)", "gllvm mem(MB)",
            "wb vs tru", "gl vs tru", "wb vs gl", "conv?"))
for (Q in c(20, 30, 50)) {
  K <- 2; J <- 30; Nj <- 15
  dat <- simulate_pfa_data(Q = Q, K = K, J = J, N_per_group = Nj, seed = Q)
  
  b_wb <- bench::mark(
    fw <<- fit_pfa_woodbury(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                            max_iter = 30, verbose = FALSE, exact_Shat = FALSE),
    iterations = 1, check = FALSE, memory = TRUE)
  t_wb <- as.numeric(b_wb$median)
  m_wb <- as.numeric(b_wb$mem_alloc) / 1024^2
  
  b_gl <- tryCatch(
    bench::mark(fr <<- fit_gllvm_ref(dat$Y, dat$X, dat$group, K = K),
                iterations = 1, check = FALSE, memory = TRUE),
    error = function(e) NULL)
  ok <- !is.null(b_gl) && fr$ok
  t_gl <- if (ok) as.numeric(b_gl$median) else NA
  m_gl <- if (ok) as.numeric(b_gl$mem_alloc) / 1024^2 else NA
  d_wb <- subspace_dist(fw$B, dat$true$B)
  d_gl <- if (ok) subspace_dist(fr$B, dat$true$B) else NA
  d_wg <- if (ok) subspace_dist(fw$B, fr$B) else NA
  
  cat(sprintf("%6d | %14.2f %14.1f %14s %14s %10.4f %10s %10s %8s\n",
              Q, t_wb, 
              m_wb,
              if (ok) sprintf("%.2f", t_gl) else "FAILED",
              if (ok) sprintf("%.1f", m_gl) else "-",
              d_wb,
              if (is.na(d_gl)) "-" else sprintf("%.4f", d_gl),
              if (is.na(d_wg)) "-" else sprintf("%.4f", d_wg),
              if (ok) ifelse(fr$converged, "yes", "NO") else "-"))
}

cat("  Q =", Q, " wb sigma2 =", round(fw$sigma2, 4),
    " true sigma2 =", round(dat$true$sigma2, 4),
    " wb vs true =", round(d_wb, 4), "\n")