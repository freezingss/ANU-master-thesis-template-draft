args <- commandArgs(trailingOnly = TRUE)
J <- as.integer(args[1])
Q <- as.integer(args[2])
seed <- as.integer(args[3])
method <- args[4]
K <- as.integer(args[5])

source("setup.R")
source("basic_functions.R")
source("sim_data.R")
source("std_multi_wbonly.R")
source("run_fit.R")

dir.create("results", showWarnings = FALSE)

result <- run_single_fit(J = J, Q = Q, K = K, seed = seed, method = method,
                         N_per_group_range = c(15, 16), M_range = c(15, 50),
                         max_iter = 80, tol = 1e-4, trace = TRUE)

fname <- sprintf("results/J%d_Q%d_seed%d_%s.rds", J, Q, seed, method)
saveRDS(result, fname)
