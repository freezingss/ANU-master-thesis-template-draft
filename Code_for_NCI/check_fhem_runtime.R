# Before dealing with Woodbury for fhem function, we check the runtime for Stage1 Estep and Mstep to decide if we should use Woodbury 
# Warm_start worth trying for fhem as well
source("fhem.R")

# fhem_backsolve_profiled <- function(lambda_hat, S_hat, tau2, ridge = 1e-3) {
#   Q <- nrow(lambda_hat)
#   J <- ncol(lambda_hat)
#   P <- diag(Q) / tau2
#   lt <- vector("list", J)
#   Fl <- vector("list", J)
#   Fl_inv <- vector("list", J)
#   nbad <- 0L
#   t_solve <- 0
#   for (j in 1:J) {
#     Sh <- S_hat[[j]]
#     t0 <- Sys.time()
#     res <- tryCatch({
#       Sh_inv <- solve(Sh)
#       Fj <- Sh_inv - P
#       Fj_reg <- Fj + diag(ridge, Q)
#       list(lt = solve(Fj_reg, Sh_inv %*% lambda_hat[, j]), Fl = Fj_reg, Fl_inv = solve(Fj_reg))
#     }, error = function(e) NULL)
#     t_solve <- t_solve + as.numeric(Sys.time() - t0, units = "secs")
#     if (is.null(res)) {
#       nbad <- nbad + 1L
#       Fl[[j]] <- diag(ridge, Q)
#       Fl_inv[[j]] <- diag(Q) / ridge
#       lt[[j]] <- rep(0, Q)
#     } else {
#       lt[[j]] <- res$lt
#       Fl[[j]] <- res$Fl
#       Fl_inv[[j]] <- res$Fl_inv
#     }
#   }
#   list(lt = lt, Fl = Fl, Fl_inv = Fl_inv, n_failed = nbad, t_solve = t_solve)
# }
# 
# fhem_gaussian_em_profiled <- function(lt, Fl, Fl_inv, Sigma0, K, max_iter = 300, tol = 1e-6) {
#   J <- length(lt)
#   q <- nrow(Sigma0)
#   Sigma <- Sigma0
#   ll <- fhem_loglik(Sigma, lt, Fl_inv)
#   trace <- ll
#   fit <- ppca_closed(Sigma, K)
#   it <- 0
#   t_sinv <- 0
#   t_estep <- 0
#   t_mstep <- 0
#   for (it in 1:max_iter) {
#     t0 <- Sys.time()
#     Sinv <- tryCatch(solve(Sigma), error = function(e) diag(1 / pmax(diag(Sigma), 1e-8)))
#     t_sinv <- t_sinv + as.numeric(Sys.time() - t0, units = "secs")
#     
#     t0 <- Sys.time()
#     S <- matrix(0, q, q)
#     for (j in 1:J) {
#       V <- solve(Fl[[j]] + Sinv)
#       # Have to be replaced by block elimination in wb_solve_exact
#       m <- V %*% (Fl[[j]] %*% lt[[j]])
#       S <- S + tcrossprod(m) + V
#     }
#     S <- S / J
#     t_estep <- t_estep + as.numeric(Sys.time() - t0, units = "secs")
#     
#     t0 <- Sys.time()
#     fit <- ppca_closed(S, K)
#     t_mstep <- t_mstep + as.numeric(Sys.time() - t0, units = "secs")
#     
#     Sigma <- fit$Sigma
#     ll_new <- fhem_loglik(Sigma, lt, Fl_inv)
#     trace <- c(trace, ll_new)
#     if (abs(ll_new - ll) < tol) break
#     ll <- ll_new
#   }
#   list(fit = fit, trace = trace, iters = it, monotone = all(diff(trace) > -1e-6),
#        t_sinv = t_sinv, t_estep = t_estep, t_mstep = t_mstep)
# }
# 
# profile_fhem <- function(K, estep_iso_fn, mstep_phi_fn, apply_PLT_fn = identity,
#                          mu0, phi0, B0, sigma2_0, ridge = 1e-3, n_outer = 3,
#                          fhem_max_iter = 300, fhem_tol = 1e-6) {
#   Q <- length(mu0)
#   mu <- mu0
#   phi <- phi0
#   B <- B0
#   sigma2 <- sigma2_0
#   tau2 <- sum(diag(tcrossprod(B))) / Q + sigma2
#   
#   timing <- data.frame(outer = integer(0), t_stage1 = numeric(0), t_mstep_phi = numeric(0),
#                        t_backsolve = numeric(0), t_sinv = numeric(0),
#                        t_fhem_estep = numeric(0), t_fhem_mstep = numeric(0),
#                        inner_iters = integer(0))
#   
#   for (outer in 1:n_outer) {
#     t0 <- Sys.time()
#     es <- estep_iso_fn(mu, phi, tau2)
#     t_stage1 <- as.numeric(Sys.time() - t0, units = "secs")
#     
#     t0 <- Sys.time()
#     mp <- mstep_phi_fn(es$lambda_hat, mu, phi)
#     mu <- mp$mu
#     phi <- mp$phi
#     t_mstep_phi <- as.numeric(Sys.time() - t0, units = "secs")
#     
#     bs <- fhem_backsolve_profiled(es$lambda_hat, es$S_hat, tau2, ridge = ridge)
#     em <- fhem_gaussian_em_profiled(bs$lt, bs$Fl, bs$Fl_inv, diag(Q) * tau2, K,
#                                     max_iter = fhem_max_iter, tol = fhem_tol)
#     
#     B <- apply_PLT_fn(em$fit$B)
#     sigma2 <- em$fit$sigma2
#     tau2 <- sum(diag(tcrossprod(B))) / Q + sigma2
#     
#     timing <- rbind(timing, data.frame(outer = outer, t_stage1 = t_stage1, t_mstep_phi = t_mstep_phi,
#                                        t_backsolve = bs$t_solve, t_sinv = em$t_sinv,
#                                        t_fhem_estep = em$t_estep, t_fhem_mstep = em$t_mstep,
#                                        inner_iters = em$iters))
#   }
#   timing
# }
# 
sim_diag <- simulate_pfa_data(Q = 400, K = 2, J = 200,
                              N_per_group_range = c(15, 16), P = 3,
                              sigma2 = 0.3, M_range = c(100, 500), seed = 1)
# 
# 
# Q <- ncol(sim_diag$Y)
# J <- max(sim_diag$group)
# th <- init_theta(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2, sigma2_init = 0.3)
# 
# estep_iso_fn <- function(mu, phi, tau2)
#   estep_wbonly(J, sim_diag$group, sim_diag$Y, sim_diag$X, sim_diag$M, mu, phi, matrix(0, Q, 2), tau2, 100, 1e-3)
# mstep_phi_fn <- function(lambda_hat, mu, phi)
#   mstep_phi_wbonly(sim_diag$Y, sim_diag$X, sim_diag$group, lambda_hat, mu, phi, lambda_phi = 0)
# 
# timing <- profile_fhem(K = 2, estep_iso_fn = estep_iso_fn, mstep_phi_fn = mstep_phi_fn,
#                        apply_PLT_fn = apply_PLT, mu0 = th$mu, phi0 = th$phi, B0 = th$B, sigma2_0 = th$sigma2,
#                        n_outer = 3, fhem_max_iter = 300)
# print(timing)
# colSums(timing[, c("t_stage1", "t_mstep_phi", "t_backsolve", "t_sinv", "t_fhem_estep", "t_fhem_mstep")])
# 
# # > print(timing) with Q = 200 J = 100 Mrange = c(50, 350)
# # outer t_stage1 t_mstep_phi t_backsolve     t_sinv t_fhem_estep t_fhem_mstep inner_iters
# # 1     1 1.271269   0.8643999   1.0430562 0.09456348    10.853158    0.3390379          23
# # 2     2 1.118330   0.6350131   1.0629511 0.06494045     7.326379    0.2331774          16
# # 3     3 0.891716   0.5937412   0.9962108 0.06343579     7.581575    0.2297750          16
# # > colSums(timing[, c("t_stage1", "t_mstep_phi", "t_backsolve", "t_sinv", "t_fhem_estep", "t_fhem_mstep")])
# # t_stage1  t_mstep_phi  t_backsolve       t_sinv t_fhem_estep t_fhem_mstep 
# # 3.2813151    2.0931542    3.1022182    0.2229397   25.7611127    0.8019903 

fit_old <- fit_pfa_fh_em_pure(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2, max_iter = 3, verbose = FALSE)
fit_new <- fit_pfa_fh_em_pure_wb(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2, max_iter = 3, verbose = FALSE)

# Check if the model setup is correct: diff should be very small
max(abs(fit_old$B - fit_new$B))
abs(fit_old$sigma2 - fit_new$sigma2)
max(abs(fit_old$trace$fh_ll - fit_new$trace$fh_ll))

# Profiling comparison
Q <- ncol(sim_diag$Y)
J <- max(sim_diag$group)
th <- init_theta(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2, sigma2_init = 0.3)

estep_iso_fn <- function(mu, phi, tau2)
  estep_wbonly_fhem_test(J, sim_diag$group, sim_diag$Y, sim_diag$X, sim_diag$M, mu, phi, matrix(0, Q, 2), tau2, 100, 1e-3)
mstep_phi_fn <- function(lambda_hat, mu, phi)
  mstep_wbonly_fhem_test(sim_diag$Y, sim_diag$X, sim_diag$group, lambda_hat, mu, phi, lambda_phi = 0)

timing_new <- profile_fhem_wb(K = 2, estep_iso_fn = estep_iso_fn, mstep_phi_fn = mstep_phi_fn,
                              apply_PLT_fn = apply_PLT, mu0 = th$mu, phi0 = th$phi, B0 = th$B, sigma2_0 = th$sigma2,
                              n_outer = 3, fhem_max_iter = 300)
print(timing_new)
colSums(timing_new[, c("t_stage1", "t_mstep_phi", "t_backsolve", "t_fhem_estep", "t_fhem_mstep", "t_fhem_loglik")])

# Check why log-likelihood differs on 1e-4
sim_small <- simulate_pfa_data(Q = 20, K = 2, J = 10,
                               N_per_group_range = c(15, 16), P = 3,
                               sigma2 = 0.3, M_range = c(10, 50), seed = 1)

fit_old3 <- fit_pfa_fh_em_pure(sim_small$Y, sim_small$X, sim_small$group, K = 2, max_iter = 1, verbose = FALSE)
fit_new3 <- fit_pfa_fh_em_pure_wb(sim_small$Y, sim_small$X, sim_small$group, K = 2, max_iter = 1, verbose = FALSE)

fit_old3$trace$fh_ll - fit_new3$trace$fh_ll
max(abs(fit_old3$B - fit_new3$B))
abs(fit_old3$sigma2 - fit_new3$sigma2)

Q <- ncol(sim_small$Y); J <- max(sim_small$group)
th <- init_theta(sim_small$Y, sim_small$X, sim_small$group, K = 2, sigma2_init = 0.3)
tau2 <- sum(diag(tcrossprod(th$B))) / Q + th$sigma2

es <- estep_wbonly_fhem_test(J, sim_small$group, sim_small$Y, sim_small$X, sim_small$M,
                             th$mu, th$phi, matrix(0, Q, 2), tau2, 100, 1e-3)
bs <- fhem_backsolve_wb(es$lambda_hat, es$d_mat, es$W_list, tau2, ridge = 1e-3)

B_test <- fit_new3$B
sigma2_test <- fit_new3$sigma2

j <- 1
W <- bs$W_list[[j]]
dt_F <- bs$d_list[[j]]
Wl <- as.numeric(W %*% bs$lt[[j]])
b <- dt_F * bs$lt[[j]] - crossprod(W, Wl)
vm <- woodbury_V_m(dt_F, W, B_test, sigma2_test, b)
quad_new <- sum(b * bs$lt[[j]]) - sum(b * vm$m)

BtB <- crossprod(B_test)
M_K <- diag(2) + BtB / sigma2_test
log_det_Sigma <- Q * log(sigma2_test) + as.numeric(determinant(M_K, logarithm = TRUE)$modulus)
logdet_new <- log_det_Sigma - bs$logdet_Freg[j] - vm$log_det_V

Fj_reg_dense <- diag(dt_F) - crossprod(W)
Fj_inv_dense <- solve(Fj_reg_dense)
Sigma_dense <- tcrossprod(B_test) + sigma2_test * diag(Q)
V_dense <- Sigma_dense + Fj_inv_dense
logdet_old <- as.numeric(determinant(V_dense, logarithm = TRUE)$modulus)
quad_old <- drop(crossprod(bs$lt[[j]], solve(V_dense, bs$lt[[j]])))

logdet_new - logdet_old
quad_new - quad_old

fit_old_r1 <- fit_pfa_fh_em_pure(sim_small$Y, sim_small$X, sim_small$group, K = 2, max_iter = 1, ridge = 1, verbose = FALSE)
fit_new_r1 <- fit_pfa_fh_em_pure_wb(sim_small$Y, sim_small$X, sim_small$group, K = 2, max_iter = 1, ridge = 1, verbose = FALSE)
fit_old_r1$trace$fh_ll - fit_new_r1$trace$fh_ll

Q <- ncol(sim_small$Y); J <- max(sim_small$group)
th <- init_theta(sim_small$Y, sim_small$X, sim_small$group, K = 2, sigma2_init = 0.3)
tau2 <- sum(diag(tcrossprod(th$B))) / Q + th$sigma2

es <- estep_wbonly_fhem_test(J, sim_small$group, sim_small$Y, sim_small$X, sim_small$M,
                             th$mu, th$phi, matrix(0, Q, 2), tau2, 100, 1e-3)
bs <- fhem_backsolve_wb(es$lambda_hat, es$d_mat, es$W_list, tau2, ridge = 1e-3)

B_test <- fit_new3$B
sigma2_test <- fit_new3$sigma2
BtB <- crossprod(B_test)
M_K <- diag(2) + BtB / sigma2_test
log_det_Sigma <- Q * log(sigma2_test) + as.numeric(determinant(M_K, logarithm = TRUE)$modulus)
Sigma_dense <- tcrossprod(B_test) + sigma2_test * diag(Q)

per_group <- data.frame(j = integer(0), Mj_total = numeric(0), min_dt = numeric(0),
                        logdet_diff = numeric(0), quad_diff = numeric(0))

for (j in 1:J) {
  W <- bs$W_list[[j]]
  dt_F <- bs$d_list[[j]]
  Wl <- as.numeric(W %*% bs$lt[[j]])
  b <- dt_F * bs$lt[[j]] - crossprod(W, Wl)
  vm <- woodbury_V_m(dt_F, W, B_test, sigma2_test, b)
  quad_new <- sum(b * bs$lt[[j]]) - sum(b * vm$m)
  logdet_new <- log_det_Sigma - bs$logdet_Freg[j] - vm$log_det_V
  
  Fj_reg_dense <- diag(dt_F) - crossprod(W)
  Fj_inv_dense <- solve(Fj_reg_dense)
  V_dense <- Sigma_dense + Fj_inv_dense
  logdet_old <- as.numeric(determinant(V_dense, logarithm = TRUE)$modulus)
  quad_old <- drop(crossprod(bs$lt[[j]], solve(V_dense, bs$lt[[j]])))
  
  per_group <- rbind(per_group, data.frame(
    j = j, Mj_total = sum(sim_small$M[sim_small$group == j]),
    min_dt = min(dt_F),
    logdet_diff = logdet_new - logdet_old,
    quad_diff = quad_new - quad_old
  ))
}
print(per_group)

Q <- ncol(sim_small$Y); J <- max(sim_small$group)
th <- init_theta(sim_small$Y, sim_small$X, sim_small$group, K = 2, sigma2_init = 0.3)
tau2 <- sum(diag(tcrossprod(th$B))) / Q + th$sigma2
es <- estep_wbonly_fhem_test(J, sim_small$group, sim_small$Y, sim_small$X, sim_small$M,
                             th$mu, th$phi, matrix(0, Q, 2), tau2, 100, 1e-3)
bs <- fhem_backsolve_wb(es$lambda_hat, es$d_mat, es$W_list, tau2, ridge = 1e-3)

ll_at_old <- fhem_loglik_wb(fit_old3$B, fit_old3$sigma2, bs$lt, bs$d_list, bs$W_list, bs$logdet_Freg)
ll_at_new <- fhem_loglik_wb(fit_new3$B, fit_new3$sigma2, bs$lt, bs$d_list, bs$W_list, bs$logdet_Freg)
ll_at_old - ll_at_new

fit_old_tight <- fit_pfa_fh_em_pure(sim_small$Y, sim_small$X, sim_small$group, K = 2, max_iter = 1, fhem_tol = 1e-10, verbose = FALSE)
fit_new_tight <- fit_pfa_fh_em_pure_wb(sim_small$Y, sim_small$X, sim_small$group, K = 2, max_iter = 1, fhem_tol = 1e-10, verbose = FALSE)
fit_old_tight$trace$fh_ll - fit_new_tight$trace$fh_ll
max(abs(fit_old_tight$B - fit_new_tight$B))

es_check <- estep_wbonly_fhem_test(J, sim_small$group, sim_small$Y, sim_small$X, sim_small$M,
                                   th$mu, th$phi, matrix(0, Q, 2), tau2, 100, 1e-3)

bs_new <- fhem_backsolve_wb(es_check$lambda_hat, es_check$d_mat, es_check$W_list, tau2, ridge = 1e-3)
bs_old <- fhem_backsolve(es_check$lambda_hat, es_check$S_hat, tau2, ridge = 1e-3)

lt_diff <- sapply(1:J, function(j) max(abs(bs_old$lt[[j]] - bs_new$lt[[j]])))
print(lt_diff)

Fl_diff <- sapply(1:J, function(j) {
  W <- es_check$W_list[[j]]
  dt_reg <- es_check$d_mat[, j] + 1e-3
  Fj_reg_new_dense <- diag(dt_reg) - crossprod(W)
  max(abs(bs_old$Fl[[j]] - Fj_reg_new_dense))
})
print(Fl_diff)

Fl_inv_diff <- sapply(1:J, function(j) {
  W <- es_check$W_list[[j]]
  dt_reg <- es_check$d_mat[, j] + 1e-3
  Fj_reg_dense <- diag(dt_reg) - crossprod(W)
  Fj_inv_fresh <- solve(Fj_reg_dense)
  max(abs(bs_old$Fl_inv[[j]] - Fj_inv_fresh))
})
print(Fl_inv_diff)

kappa_vals <- sapply(1:J, function(j) {
  W <- es_check$W_list[[j]]
  dt_reg <- es_check$d_mat[, j] + 1e-3
  kappa(diag(dt_reg) - crossprod(W))
})
print(kappa_vals)

cbind(Fl_inv_diff, kappa_vals)