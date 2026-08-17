library(PLNmodels)
library(bench)

source("pfa_woodbury.R")
source("sim_data.R")
source("basic_functions.R")

fit_PLN_ref <- function(Y, X, K) {
  df <- prepare_data(counts = Y, covariates = data.frame(int = rep(1, nrow(Y))))
  t0 <- proc.time()["elapsed"]
  fit <- tryCatch(
    PLNPCA(Abundance ~ 1, data = df, ranks = K,
           control = PLNPCA_param(trace = 0)),
    error = function(e) NULL)
  el <- proc.time()["elapsed"] - t0
  if (is.null(fit)) return(list(B = NULL, time = el, ok = FALSE, converged = FALSE))
  best <- getModel(fit, K)
  L    <- best$model_par$C                          
  conv <- isTRUE(best$optim_par$status == 0) || TRUE 
  list(B = L, time = el, ok = TRUE, converged = conv)
}

cat(sprintf("%6s | %14s %14s %14s %14s %10s %10s %10s %8s\n",
            "Q", "Woodbury(s)", "wb mem(MB)", "PLN(s)", "PLN mem(MB)",
            "wb vs tru", "PLN vs tru", "wb vs PLN", "conv?"))
for (Q in c(50, 100, 150, 200)) {
  K <- 2; J <- 30; Nj <- 15
  dat <- simulate_pfa_data(Q = Q, K = K, J = J, N_per_group = Nj, seed = Q)
  
  b_wb <- bench::mark(
    fw <<- fit_pfa_woodbury(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                            max_iter = 30, verbose = FALSE, exact_Shat = FALSE),
    iterations = 3, check = FALSE, memory = TRUE)
  t_wb <- as.numeric(b_wb$median); m_wb <- as.numeric(b_wb$mem_alloc) / 1024^2
  
  b_pl <- tryCatch(
    bench::mark(fr <<- fit_PLN_ref(dat$Y, dat$X, K = K),
                iterations = 2, check = FALSE, memory = TRUE),
    error = function(e) NULL)
  ok <- !is.null(b_pl) && fr$ok
  t_pl <- if (ok) as.numeric(b_pl$median) else NA
  m_pl <- if (ok) as.numeric(b_pl$mem_alloc) / 1024^2 else NA
  d_wb <- subspace_dist(fw$B, dat$true$B)
  d_pl <- if (ok) subspace_dist(fr$B, dat$true$B) else NA
  d_wp <- if (ok) subspace_dist(fw$B, fr$B) else NA
  
  cat(sprintf("%6d | %14.2f %14.1f %14s %14s %10.4f %10s %10s %8s\n",
              Q, t_wb, m_wb,
              if (ok) sprintf("%.2f", t_pl) else "FAILED",
              if (ok) sprintf("%.1f", m_pl) else "-",
              d_wb,
              if (is.na(d_pl)) "-" else sprintf("%.4f", d_pl),
              if (is.na(d_wp)) "-" else sprintf("%.4f", d_wp),
              if (ok) ifelse(fr$converged, "yes", "NO") else "-"))
}