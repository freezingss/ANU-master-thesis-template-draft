source("basic_functions.R")
source("pfa_base_dense.R")
source("pfa_woodbury_only.R")
source("pfa_woodbury_lam_corrected.R")
source("pfa_debias.R")
source("pfa_fhem.R")
source("pfa_fhem_woodbury.R")
source("sim_data.R")

# res <- run_grid(Js = c(25, 50, 75), 
#                 Qs = c(25, 50, 75), 
#                 seeds = 1, 
#                 K = 2, 
#                 Nj = 15,
#                 sigma2_true = 0.3, 
#                 M_range = c(50, 150), 
#                 max_iter = 80, 
#                 tol = 1e-4,
#                 methods = c("std", "wb", "corrected", "decorr", "debias", "fh_em", "fh_em_pure", "glmmTMB"))
# summarize_grid(res)

# Version 2: run the test with correction Pi and output Pi*Sigma*Pi
deflate_cols <- function(B) {
  if (is.null(B)) return(NULL)
  B - matrix(colMeans(B), nrow(B), ncol(B), byrow = TRUE)
}

u_frac <- function(B) {
  if (is.null(B)) return(NA_real_)
  Q <- nrow(B)
  sqrt(sum(colSums(B)^2) / Q) / sqrt(sum(B^2))
}

# GLF Results Table: compared with isotropic, existing packages and anisotropic methods, given discretely different phenomenons

res <- run_grid(Js = c(50, 75, 100),
                Qs = c(75, 100, 125, 150),
                seeds = 1,
                K = 2,
                Nj = 15,
                sigma2_true = 0.3,
                M_range = c(75, 150),
                max_iter = 80,
                tol = 1e-4,
                methods = c("wb", "corrected", "fhem", "fhem_wb", "glmmTMB"))
# summarize_grid(res)
saveRDS(res, file = "res_grid_14092026v3.rds")

# Run the code with improved glmmTMB package rr()
res_offset <- run_grid(Js = c(50),
                Qs = c(75),
                seeds = 1,
                K = 2,
                Nj = 15,
                sigma2_true = 0.3,
                M_range = c(75, 150),
                max_iter = 80,
                tol = 1e-4,
                methods = c("wb", "corrected", "fhem", "fhem_wb", "glmmTMB"))

# Check 1
res <- readRDS("res_grid_14092026v1.rds")
aggregate(u_frac ~ method, data = res, FUN = function(x) c(mean = mean(x), max = max(x)))
aggregate(u_frac ~ method + J + Q, data = res, FUN = mean)
# Problem happens, the GLF value of wb & correction not numerically 0, implies that our projection doesn't work well
# Try: deflate the columns of estimated B and then do the PPCA to avoid accumulating the error


# Check 2
res_check2 <- run_grid(Js = c(100),
                Qs = c(150),
                seeds = 1,
                K = 2,
                Nj = 15,
                sigma2_true = 0.3,
                M_range = c(75, 150),
                max_iter = 80,
                tol = 1e-4,
                methods = c("wb", "corrected", "fhem", "fhem_wb", "glmmTMB"))
# summarize_grid(res)
saveRDS(res_check2, file = "res_grid_14092026v2.rds")

Bs <- attr(res, "B_store")
energy_check <- do.call(rbind, lapply(Bs, function(rec) {
  data.frame(J = rec$J, Q = rec$Q, seed = rec$seed, method = rec$method,
             energy_hat = sum(rec$B_hat_deflated^2),
             energy_true = sum(rec$B_true_deflated^2),
             ratio = sum(rec$B_hat_deflated^2) / sum(rec$B_true_deflated^2))
}))
energy_check[order(energy_check$J, energy_check$Q, energy_check$method), ]

# Check 3: new ppca
