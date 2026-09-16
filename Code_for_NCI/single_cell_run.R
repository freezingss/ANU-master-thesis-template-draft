# Single grid for different seeds
run_single_fit <- function(J, Q, seed, method, K,
                           P = 3, N_per_group_range = c(10, 20),
                           M_range = c(200, 500), sigma2_true = 0.3) {
  dat <- simulate_pfa_data(Q = Q, K = K, J = J, N_per_group_range = N_per_group_range,
                           P = P, sigma2 = sigma2_true, M_range = M_range, seed = seed)
  
  fit <- switch(method,
                wb = fit_pfa_wbonly(dat$Y, dat$X, dat$group, K,
                                    B_true = dat$true$B, trace = TRUE),
                corrected = fit_pfa_woodbury_lam_corr(dat$Y, dat$X, dat$group, K,
                                                      B_true = dat$true$B, use_lambda_correction = TRUE, trace = TRUE),
                fhem = fit_pfa_fhem(dat$Y, dat$X, dat$group, K,
                                    B_true = dat$true$B, trace = TRUE),
                glmmTMB = fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K),
                stop(sprintf("unknown method: %s", method)))
  
  list(config = list(J = J, Q = Q, seed = seed, method = method, K = K),
       truth = dat$true, fit = fit)
}

args <- commandArgs(trailingOnly = TRUE)
J <- as.integer(args[1])
Q <- as.integer(args[2])
seed <- as.integer(args[3])
method <- args[4]
K <- as.integer(args[5])

source("setup.R")
source("basic_functions.R")
source("sim_data.R")
source("std_multinomial_wb.R")
source("std_multinomial_wb_lam_corrected.R")
source("fhem_wb.R")
source("existing_packages.R")

result <- run_single_fit(J, Q, seed, method, K)
saveRDS(result, sprintf("results/J%d_Q%d_seed%d_%s.rds", J, Q, seed, method))