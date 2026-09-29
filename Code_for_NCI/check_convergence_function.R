# Check the convergence issue for std_multinomial_wb
source("setup.R")
source("basic_functions.R")
source("std_multinomial_wb.R")
source("run_and_plot.R") # run_grid function defined here
source("fit_pfa_wbonly_warm_start.R")

# Check some criterion

res_grid_multi_wb_check <- run_grid(Js = c(50),
                     Qs = c(50),
                     seeds = 1,
                     K = 2, 
                     Nj = 15, 
                     sigma2_true = 0.3,
                     M_range = c(75, 150),
                     max_iter = 80, 
                     tol = 1e-4,
                     methods = c("wb")
                     )

# res_repsmulti_wb_check <- run_grid(Js = c(200),
#                      Qs = c(400),
#                      seeds = 2:4,
#                      K = 2, 
#                      Nj = 15, 
#                      sigma2_true = 0.3,
#                      M_range = c(200, 600),
#                      max_iter = 80, 
#                      tol = 1e-4,
#                      methods = c("wb")
#                      )

# 1. Monotonicity: theoretically non-monotone because of not exact EM algorithm, introducing approximation

check_monotonicity <- function(Js, Qs, seed, K, Nj, sigma2_true, M_range, max_iter, tol, P = 3) {
  out <- data.frame()
  for (J in Js) {
    for (Q in Qs) {
      sim <- simulate_pfa_data(Q = Q, 
                               K = K, 
                               J = J,
                               N_per_group_range = c(Nj, Nj), 
                               P = P,
                               sigma2 = sigma2_true, 
                               M_range = M_range, 
                               seed = seed)
      fit <- fit_pfa_wbonly(sim$Y, 
                            sim$X, 
                            sim$group, K,
                            max_iter = max_iter, 
                            tol = tol,
                            trace = TRUE, 
                            B_true = sim$true$B)
      peak_iter <- which.max(fit$log_evidence)
      bdist_min_iter <- which.min(fit$B_dist_trace)
      out <- rbind(out, data.frame(
        J = J, Q = Q,
        n_iter = fit$iterations,
        stop_reason = fit$stop_reason,
        is_monotone = fit$is_monotone,
        first_decrease_iter = fit$first_decrease_iter,
        peak_ev_iter = peak_iter,
        bdist_min_iter = bdist_min_iter,
        bdist_final = tail(fit$B_dist_trace, 1),
        bdist_at_peak = fit$B_dist_trace[peak_iter]
      ))
    }
  }
  out
}

mono_check <- check_monotonicity(Js = c(50), 
                                 Qs = c(50), 
                                 seed = 1, 
                                 K = 2, 
                                 Nj = 15,
                                 sigma2_true = 0.3, 
                                 M_range = c(25, 75),
                                 max_iter = 200, 
                                 tol = 1e-4
                                 )
print(mono_check)
cat(sprintf("bdist_at_peak < bdist_final: %d / %d cell\n",
            sum(mono_check$bdist_at_peak < mono_check$bdist_final), nrow(mono_check)))

# 2. Decrease shape

sim75 <- simulate_pfa_data(Q = 50, 
                           K = 2, 
                           J = 50,
                           N_per_group_range = c(15, 15), 
                           P = 3,
                           sigma2 = 0.3, 
                           M_range = c(25, 75), 
                           seed = 1)
fit_long <- fit_pfa_wbonly(sim75$Y, 
                           sim75$X, 
                           sim75$group, 
                           K = 2,
                           max_iter = 400, 
                           tol = 1e-4,
                           trace = TRUE, 
                           B_true = sim75$true$B)

fit_long$converged
fit_long$iterations
d <- diff(fit_long$log_evidence)
tail(d, 20)

# Estimated: how many steps reduce the absolute value of log-likelihood lower than the tolerance?

estimate_geometric_convergence <- function(log_ev, tol = 1e-4, window = 30) {
  d <- diff(log_ev)
  n_d <- length(d)
  idx <- max(1, n_d - window + 1):n_d
  fit <- lm(log(d[idx]) ~ idx)
  r <- exp(coef(fit)[2])
  d_last <- d[n_d]
  remaining_gap <- d_last * r / (1 - r)
  iters_needed <- if (d_last > tol) ceiling(log(tol / d_last) / log(r)) else 0
  list(r = as.numeric(r), d_last = d_last,
       remaining_gap = remaining_gap,
       extrapolated_log_ev = log_ev[length(log_ev)] + remaining_gap,
       iters_to_tol = iters_needed)
}

estimate_geometric_convergence(fit_long$log_evidence, tol = 1e-4, window = 30)

# Considering to use relative absolute stopping criterion, since a small test with max_iter = 400
# keeps unconverged.

tail(fit_long$sigma2_trace, 20)
diff(tail(fit_long$sigma2_trace, 21))
tail(fit_long$B_dist_trace, 20)
diff(tail(fit_long$B_dist, 21))


# Held-out try

fit_pfa_wbonly_heldout <- function(Y, X, group, K,
                                   M = rowSums(Y),
                                   heldout_frac = 0.2,
                                   max_iter = 60,
                                   lambda_phi = 0, sigma2_init = 0.3,
                                   estep_max_iter = 100, estep_gtol = 1e-3,
                                   patience = 15, seed_split = 1,
                                   trace = FALSE, B_true = NULL) {
  
  J <- max(group)
  set.seed(seed_split)
  heldout_groups <- sample(1:J, size = round(J * heldout_frac))
  train_groups <- setdiff(1:J, heldout_groups)
  
  train_idx <- which(group %in% train_groups)
  heldout_idx <- which(group %in% heldout_groups)
  
  Y_tr <- Y[train_idx, , drop = FALSE]
  X_tr <- X[train_idx, , drop = FALSE]
  group_tr <- match(group[train_idx], train_groups)
  M_tr <- M[train_idx]
  J_tr <- length(train_groups)
  
  Y_ho <- Y[heldout_idx, , drop = FALSE]
  X_ho <- X[heldout_idx, , drop = FALSE]
  group_ho <- match(group[heldout_idx], heldout_groups)
  M_ho <- M[heldout_idx]
  J_ho <- length(heldout_groups)
  
  Q <- ncol(Y)
  P <- ncol(X)
  stopifnot(all(X[, 1] == 1))
  
  avg_prop <- colMeans(Y_tr / pmax(rowSums(Y_tr), 1))
  mu <- log(avg_prop + 1e-8)
  mu <- mu - mean(mu)
  phi <- matrix(0, P, Q)
  
  gm <- matrix(0, J_tr, Q)
  for (j in 1:J_tr) {
    idx <- which(group_tr == j)
    if (length(idx) > 0) {
      gs <- colSums(Y_tr[idx, , drop = FALSE])
      gm[j, ] <- log((gs + 1e-5) / (sum(gs) + Q * 1e-5))
    }
  }
  gm_c <- sweep(gm, 2, colMeans(gm))
  sv0 <- svd(t(gm_c), nu = K, nv = K)
  B <- apply_PLT(sv0$u %*% diag(pmax(sv0$d[1:K] * 0.5, 0.1), K))
  sigma2 <- sigma2_init
  
  log_ev_train <- numeric(max_iter)
  log_ev_heldout <- numeric(max_iter)
  if (trace && !is.null(B_true)) B_dist_trace <- numeric(max_iter)
  
  best_ho <- -Inf
  best_iter <- NA_integer_
  best_state <- NULL
  no_improve <- 0
  em <- 0
  
  for (em in 1:max_iter) {
    
    es_tr <- estep_wbonly(J_tr, group_tr, Y_tr, X_tr, M_tr, mu, phi, B, sigma2,
                          max_iter = estep_max_iter, estep_gtol = estep_gtol)
    log_ev_train[em] <- es_tr$lp_total - 0.5 * J_tr * es_tr$log_det_Sigma + 0.5 * es_tr$ld_S_total
    
    es_ho <- estep_wbonly(J_ho, group_ho, Y_ho, X_ho, M_ho, mu, phi, B, sigma2,
                          max_iter = estep_max_iter, estep_gtol = estep_gtol)
    log_ev_heldout[em] <- es_ho$lp_total - 0.5 * J_ho * es_ho$log_det_Sigma + 0.5 * es_ho$ld_S_total
    
    if (log_ev_heldout[em] > best_ho) {
      best_ho <- log_ev_heldout[em]
      best_iter <- em
      best_state <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2)
      no_improve <- 0
    } else {
      no_improve <- no_improve + 1
    }
    
    mp <- mstep_phi_wbonly(Y_tr, X_tr, group_tr, es_tr$lambda_hat, mu, phi, lambda_phi = lambda_phi)
    mu <- mp$mu
    phi <- mp$phi
    
    S_obs <- tcrossprod(es_tr$lambda_hat) / J_tr
    for (j in 1:J_tr) S_obs <- S_obs + es_tr$S_hat[[j]] / J_tr
    rt <- ppca_closed(S_obs, K)
    B <- apply_PLT(rt$B)
    sigma2 <- rt$sigma2
    
    if (trace && !is.null(B_true)) B_dist_trace[em] <- subspace_distance(B, B_true)
    
    if (no_improve >= patience) break
  }
  
  em_final <- em
  out <- list(mu = best_state$mu, phi = best_state$phi,
              B = best_state$B, sigma2 = best_state$sigma2,
              log_ev_train = log_ev_train[1:em_final],
              log_ev_heldout = log_ev_heldout[1:em_final],
              best_iter = best_iter, iterations = em_final)
  if (trace && !is.null(B_true)) out$B_dist_trace <- B_dist_trace[1:em_final]
  out
}

fit_cv <- fit_pfa_wbonly_heldout(sim75$Y, sim75$X, sim75$group, K = 2,
                                 heldout_frac = 0.2, max_iter = 200, patience = 15,
                                 trace = TRUE, B_true = sim75$true$B)

fit_cv$iterations
fit_cv$best_iter
which.min(fit_cv$B_dist_trace)
subspace_distance(fit_cv$B, sim75$true$B)

# Correct the stable converged point try

fit_full <- fit_pfa_wbonly(sim75$Y, sim75$X, sim75$group, K = 2,
                           max_iter = 200, tol = 1e-4,
                           trace = TRUE, B_true = sim75$true$B)

J <- max(sim75$group)
Q <- ncol(sim75$Y)

cc <- apply_lambda_correction(Q = Q, J = J,
                              lambda_hat = fit_full$lambda_hat, S_hat = fit_full$S_hat,
                              Y = sim75$Y, X = sim75$X, group = sim75$group,
                              M = rowSums(sim75$Y), mu = fit_full$mu, phi = fit_full$phi,
                              corr_max_rel = 1)

S_obs_corr <- tcrossprod(cc$lambda_corrected) / J

for (j in 1:J) S_obs_corr <- S_obs_corr + fit_full$S_hat[[j]] / J
rt_corr <- ppca_closed(S_obs_corr, 2)
B_corr <- apply_PLT(rt_corr$B)

subspace_distance(fit_full$B, sim75$true$B)
subspace_distance(B_corr, sim75$true$B)

floor_val <- bbp_floor(sim75$true$B, sim75$true$sigma2, J = 50)
floor_val
subspace_distance(fit_full$B, sim75$true$B)
subspace_distance(B_corr, sim75$true$B)

S_true_lambda <- tcrossprod(sim75$true$lambda) / J
rt_true_lambda <- ppca_closed(S_true_lambda, 2)
B_true_pca <- apply_PLT(rt_true_lambda$B)
subspace_distance(B_true_pca, sim75$true$B)

# Compared true - fit - pca - bbp

seeds <- 2:5

floor_vec <- numeric(length(seeds))
pca_dist_vec <- numeric(length(seeds))
fit_dist_vec <- numeric(length(seeds))
iters_vec <- numeric(length(seeds))
stop_reason_vec <- character(length(seeds))

for (i in seq_along(seeds)) {
  seed <- seeds[i]
  
  sim75_bigM <- simulate_pfa_data(Q = 50, K = 2, J = 50,
                                  N_per_group_range = c(15, 15), P = 3,
                                  sigma2 = 0.3, M_range = c(60, 120), seed = seed)
  
  fit_bigM <- fit_pfa_wbonly(sim75_bigM$Y, sim75_bigM$X, sim75_bigM$group, K = 2,
                             max_iter = 200, tol = 1e-4,
                             trace = TRUE, B_true = sim75_bigM$true$B)
  
  S_true_lambda <- tcrossprod(sim75_bigM$true$lambda) / 50
  rt_true <- ppca_closed(S_true_lambda, 2)
  B_true_pca <- apply_PLT(rt_true$B)
  
  floor_vec[i] <- bbp_floor(sim75_bigM$true$B, sim75_bigM$true$sigma2, J = 50)
  pca_dist_vec[i] <- subspace_distance(B_true_pca, sim75_bigM$true$B)
  fit_dist_vec[i] <- subspace_distance(fit_bigM$B, sim75_bigM$true$B)
  iters_vec[i] <- fit_bigM$iterations
  stop_reason_vec[i] <- fit_bigM$stop_reason
}

summary_tab <- data.frame(seed = seeds, floor = floor_vec, pca_dist = pca_dist_vec,
                          fit_dist = fit_dist_vec, iterations = iters_vec,
                          stop_reason = stop_reason_vec)
print(summary_tab)
colMeans(summary_tab[, c("floor", "pca_dist", "fit_dist")])

fit_diag <- fit_pfa_wbonly(sim75_bigM$Y, sim75_bigM$X, sim75_bigM$group, K = 2,
                           max_iter = 200, tol = 1e-4, trace = TRUE, B_true = sim75_bigM$true$B)
tail(fit_diag$B_shift_trace, 20)
tail(fit_diag$sigma2_shift_trace, 20)
tail(fit_diag$r_hat_trace, 20)
tail(fit_diag$remaining_gap_trace, 20)

# Monte Carlo rep for lambda and S_hat

check_S_hat_calibration <- function(j_index, sim, R = 300) {
  idx <- which(sim$group == j_index)
  X_j <- sim$X[idx, , drop = FALSE]
  M_j <- sim$M[idx]
  lambda_true_j <- sim$true$lambda[, j_index]
  mu <- sim$true$mu
  phi <- sim$true$phi
  B <- sim$true$B
  sigma2 <- sim$true$sigma2
  Q <- ncol(sim$Y)
  K <- ncol(B)
  
  fixed <- matrix(rep(mu, each = length(idx)), length(idx), Q) + X_j %*% phi
  eta_true <- sweep(fixed, 2, lambda_true_j, "+")
  pi_true <- row_softmax(eta_true)
  
  BtB <- crossprod(B)
  M_K <- diag(K) + BtB / sigma2
  M_K_inv <- solve(M_K)
  log_det_MK <- 2 * sum(log(diag(chol(M_K))))
  
  lambda_hat_rep <- matrix(0, Q, R)
  S_hat_diag_rep <- matrix(0, Q, R)
  
  for (r in 1:R) {
    Y_j_rep <- t(sapply(seq_along(idx), function(i) rmultinom(1, M_j[i], pi_true[i, ])))
    res <- laplace_lambda_j_wbonly(Y_j_rep, X_j, M_j, mu, phi, B, sigma2,
                                   BtB = BtB, M_K = M_K, M_K_inv = M_K_inv,
                                   log_det_MK = log_det_MK)
    lambda_hat_rep[, r] <- res$lambda_hat
    S_hat_diag_rep[, r] <- diag(res$S_hat)
  }
  
  empirical_var <- apply(lambda_hat_rep, 1, var)
  model_var <- rowMeans(S_hat_diag_rep)
  
  list(bias = rowMeans(lambda_hat_rep) - lambda_true_j,
       empirical_var = empirical_var, model_var = model_var,
       ratio = empirical_var / model_var)
}

# Redefine sim75_bigM to make sure only M_range differs

sim75_bigM <- simulate_pfa_data(Q = 50, K = 2, J = 50,
                                N_per_group_range = c(15, 15), P = 3,
                                sigma2 = 0.3, M_range = c(60, 120), seed = 1)

# Single test
identical(sim75$true$lambda, sim75_bigM$true$lambda)

res_smallM <- check_S_hat_calibration(j_index = 1, sim = sim75, R = 30)
res_bigM <- check_S_hat_calibration(j_index = 1, sim = sim75_bigM, R = 30)

mean(res_smallM$ratio)
mean(res_bigM$ratio)

# More seeds rep

check_S_hat_calibration_multi <- function(sim, j_indices, R = 300) {
  out <- lapply(j_indices, function(j) {
    res <- check_S_hat_calibration(j_index = j, sim = sim, R = R)
    data.frame(j = j, ratio = res$ratio, bias = res$bias,
               empirical_var = res$empirical_var, model_var = res$model_var)
  })
  do.call(rbind, out)
}

cal_smallM <- check_S_hat_calibration_multi(sim75, j_indices = 1:10, R = 300)
cal_bigM <- check_S_hat_calibration_multi(sim75_bigM, j_indices = 1:10, R = 300)

mean(cal_smallM$ratio)
mean(cal_bigM$ratio)

summary(cal_smallM$ratio)
summary(cal_bigM$ratio)

# Exact test: effect from multinomial

check_S_hat_calibration <- function(j_index, sim, R = 300, M_target = NULL) {
  idx <- which(sim$group == j_index)
  X_j <- sim$X[idx, , drop = FALSE]
  M_j <- if (is.null(M_target)) sim$M[idx] else rep(M_target, length(idx))
  lambda_true_j <- sim$true$lambda[, j_index]
  mu <- sim$true$mu
  phi <- sim$true$phi
  B <- sim$true$B
  sigma2 <- sim$true$sigma2
  Q <- ncol(sim$Y)
  K <- ncol(B)
  
  fixed <- matrix(rep(mu, each = length(idx)), length(idx), Q) + X_j %*% phi
  eta_true <- sweep(fixed, 2, lambda_true_j, "+")
  pi_true <- row_softmax(eta_true)
  
  BtB <- crossprod(B)
  M_K <- diag(K) + BtB / sigma2
  M_K_inv <- solve(M_K)
  log_det_MK <- 2 * sum(log(diag(chol(M_K))))
  
  lambda_hat_rep <- matrix(0, Q, R)
  S_hat_diag_rep <- matrix(0, Q, R)
  grad_norm_rep <- numeric(R)
  n_iter_rep <- numeric(R)
  
  for (r in 1:R) {
    Y_j_rep <- t(sapply(seq_along(idx), function(i) rmultinom(1, M_j[i], pi_true[i, ])))
    res <- laplace_lambda_j_wbonly(Y_j_rep, X_j, M_j, mu, phi, B, sigma2,
                                   BtB = BtB, M_K = M_K, M_K_inv = M_K_inv,
                                   log_det_MK = log_det_MK)
    lambda_hat_rep[, r] <- res$lambda_hat
    S_hat_diag_rep[, r] <- diag(res$S_hat)
    grad_norm_rep[r] <- res$grad_norm
    n_iter_rep[r] <- res$n_iter
  }
  
  var_per_q <- apply(lambda_hat_rep, 1, var)
  mse_per_q <- rowMeans((lambda_hat_rep - lambda_true_j)^2)
  model_var <- rowMeans(S_hat_diag_rep)
  
  list(bias = rowMeans(lambda_hat_rep) - lambda_true_j,
       ratio_var = var_per_q / model_var, ratio_mse = mse_per_q / model_var,
       grad_norm = grad_norm_rep, n_iter = n_iter_rep)
}

M_grid <- c(10, 30, 60, 100, 150, 300, 600, 1200, 3000, 6000)
j_indices <- 1:10

sweep_out <- lapply(M_grid, function(Mt) {
  res <- lapply(j_indices, function(j) {
    r <- check_S_hat_calibration(j_index = j, sim = sim75, R = 300, M_target = Mt)
    data.frame(M_target = Mt, j = j, ratio_var = mean(r$ratio_var), ratio_mse = mean(r$ratio_mse),
               grad_norm_max = max(r$grad_norm), n_iter_max = max(r$n_iter))
  })
  do.call(rbind, res)
})
sweep_tab <- do.call(rbind, sweep_out)

agg <- aggregate(cbind(ratio_var, ratio_mse, grad_norm_max, n_iter_max) ~ M_target, data = sweep_tab, FUN = mean)
print(agg)

# Abnormal: try projection

check_S_hat_calibration_perp <- function(j_index, sim, R = 300, M_target = NULL) {
  idx <- which(sim$group == j_index)
  X_j <- sim$X[idx, , drop = FALSE]
  M_j <- if (is.null(M_target)) sim$M[idx] else rep(M_target, length(idx))
  lambda_true_j <- sim$true$lambda[, j_index]
  mu <- sim$true$mu
  phi <- sim$true$phi
  B <- sim$true$B
  sigma2 <- sim$true$sigma2
  Q <- ncol(sim$Y)
  K <- ncol(B)
  
  V <- qr.Q(qr(matrix(1, Q, 1)), complete = TRUE)[, 2:Q, drop = FALSE]
  lambda_true_perp <- as.numeric(t(V) %*% lambda_true_j)
  
  fixed <- matrix(rep(mu, each = length(idx)), length(idx), Q) + X_j %*% phi
  eta_true <- sweep(fixed, 2, lambda_true_j, "+")
  pi_true <- row_softmax(eta_true)
  
  BtB <- crossprod(B)
  M_K <- diag(K) + BtB / sigma2
  M_K_inv <- solve(M_K)
  log_det_MK <- 2 * sum(log(diag(chol(M_K))))
  
  lambda_hat_perp_rep <- matrix(0, Q - 1, R)
  S_hat_perp_diag_rep <- matrix(0, Q - 1, R)
  
  for (r in 1:R) {
    Y_j_rep <- t(sapply(seq_along(idx), function(i) rmultinom(1, M_j[i], pi_true[i, ])))
    res <- laplace_lambda_j_wbonly(Y_j_rep, X_j, M_j, mu, phi, B, sigma2,
                                   BtB = BtB, M_K = M_K, M_K_inv = M_K_inv,
                                   log_det_MK = log_det_MK)
    lambda_hat_perp_rep[, r] <- as.numeric(t(V) %*% res$lambda_hat)
    S_hat_perp_diag_rep[, r] <- diag(t(V) %*% res$S_hat %*% V)
  }
  
  mse_perp <- rowMeans((lambda_hat_perp_rep - lambda_true_perp)^2)
  var_perp <- apply(lambda_hat_perp_rep, 1, var)
  model_var_perp <- rowMeans(S_hat_perp_diag_rep)
  
  list(ratio_var_perp = var_perp / model_var_perp,
       ratio_mse_perp = mse_perp / model_var_perp)
}

M_grid <- c(10, 30, 60, 100, 150, 300, 600, 1200, 3000, 6000)
j_indices <- 1:10

sweep_out <- lapply(M_grid, function(Mt) {
  res <- lapply(j_indices, function(j) {
    r <- check_S_hat_calibration_perp(j_index = j, sim = sim75, R = 300, M_target = Mt)
    data.frame(M_target = Mt, j = j, ratio_var_perp = mean(r$ratio_var_perp),
               ratio_mse_perp = mean(r$ratio_mse_perp))
  })
  do.call(rbind, res)
})
sweep_tab_perp <- do.call(rbind, sweep_out)
agg_perp <- aggregate(cbind(ratio_var_perp, ratio_mse_perp) ~ M_target, data = sweep_tab_perp, FUN = mean)
print(agg_perp)

fit_diag <- fit_pfa_wbonly(sim75$Y, sim75$X, sim75$group, K = 2,
                           max_iter = 80, tol = 1e-4, trace = TRUE, B_true = sim75$true$B)
range(fit_diag$mu_bar_trace)
apply(fit_diag$phi_bar_trace, 1, range)
plot(fit_diag$mu_bar_trace, type = "l")

colSums(sim75$true$B)
colMeans(sim75$true$B)

check_monotonicity <- function(Js, Qs, seeds, K, Nj, sigma2_true, M_range, max_iter, tol,
                               P = 3, patience = 15, burn_in = 15, rate_window = 10) {
  out <- data.frame()
  for (J in Js) {
    for (Q in Qs) {
      for (seed in seeds) {
        sim <- simulate_pfa_data(Q = Q, K = K, J = J,
                                 N_per_group_range = c(Nj, Nj), P = P,
                                 sigma2 = sigma2_true, M_range = M_range, seed = seed)
        fit <- fit_pfa_wbonly(sim$Y, sim$X, sim$group, K,
                              max_iter = max_iter, tol = tol,
                              patience = patience, burn_in = burn_in, rate_window = rate_window,
                              trace = TRUE, B_true = sim$true$B)
        peak_iter <- which.max(fit$log_evidence)
        bdist_min_iter <- which.min(fit$B_dist_trace)
        out <- rbind(out, data.frame(
          J = J, Q = Q, seed = seed,
          n_iter = fit$iterations,
          stop_reason = fit$stop_reason,
          is_monotone = fit$is_monotone,
          first_decrease_iter = fit$first_decrease_iter,
          peak_ev_iter = peak_iter,
          bdist_min_iter = bdist_min_iter,
          bdist_final = tail(fit$B_dist_trace, 1),
          bdist_at_peak = fit$B_dist_trace[peak_iter]
        ))
      }
    }
  }
  out
}

mono_check2 <- check_monotonicity(Js = c(50), Qs = c(75, 100, 125, 150), seeds = 1:3,
                                  K = 2, Nj = 15, sigma2_true = 0.3, M_range = c(75, 150),
                                  max_iter = 80, tol = 1e-4)
print(mono_check2)

# Why three criterion do not make sense?

sim_diag <- simulate_pfa_data(Q = 100, K = 2, J = 50,
                              N_per_group_range = c(15, 15), P = 3,
                              sigma2 = 0.3, M_range = c(75, 150), seed = 2)
fit_diag <- fit_pfa_wbonly(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2,
                           max_iter = 80, tol = 1e-3,
                           trace = TRUE, B_true = sim_diag$true$B)

tail(fit_diag$mu_shift_trace, 20)
tail(fit_diag$phi_shift_trace, 20)
tail(fit_diag$no_stall_trace, 20)
fit_diag$stop_reason
fit_diag$iterations

# Add the automatic stopping standard for four parameters

fit_diag_warm_start <- fit_pfa_wbonly_warm_start(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2,
                                                 max_iter = 80, tol_B = 1e-4, tol_sigma2 = 1e-4,
                                                 trace = TRUE, B_true = sim_diag$true$B)

# fit_diag_warm_start$tol_mu
# fit_diag_warm_start$tol_phi

fit_diag_warm_start$tol_sigma2
fit_diag_warm_start$tol_B
fit_diag_warm_start$stop_reason
fit_diag_warm_start$iterations

# tail(fit_diag_warm_start$mu_r_hat_trace, 10)
# tail(fit_diag_warm_start$mu_remaining_trace, 10)

# tail(fit_diag_warm_start$mu_shift_trace, 20)
# tail(fit_diag_warm_start$phi_shift_trace, 20)
# tail(fit_diag_warm_start$no_stall_trace, 20)
# fit_diag_warm_start$stop_reason
# fit_diag_warm_start$iterations

# fit_gtol3 <- fit_pfa_wbonly_warm_start(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2,
#                                        max_iter = 80, tol = 1e-4, estep_gtol = 1e-3,
#                                        trace = TRUE, B_true = sim_diag$true$B)
# fit_gtol6 <- fit_pfa_wbonly_warm_start(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2,
#                                        max_iter = 80, tol = 1e-4, estep_gtol = 1e-6,
#                                        trace = TRUE, B_true = sim_diag$true$B)
# 
# tail(fit_gtol3$mu_shift_trace, 20)
# tail(fit_gtol6$mu_shift_trace, 20)

# Compare the estimation for probability from different mu and different iterations, maximal difference from two probabilities on 1e-4

compare_fitted_pi <- function(fit_a, fit_b, X, group) {
  eta_a <- sweep(X %*% fit_a$phi, 2, fit_a$mu, "+") + t(fit_a$lambda_hat[, group, drop = FALSE])
  eta_b <- sweep(X %*% fit_b$phi, 2, fit_b$mu, "+") + t(fit_b$lambda_hat[, group, drop = FALSE])
  pi_a <- row_softmax(eta_a)
  pi_b <- row_softmax(eta_b)
  max(abs(pi_a - pi_b))
}

sim_diag <- simulate_pfa_data(Q = 200, K = 2, J = 50,
                              N_per_group_range = c(15, 16), P = 3,
                              sigma2 = 0.3, M_range = c(150, 300), seed = 1)

fit_short <- fit_pfa_wbonly_warm_start(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2,
                                       max_iter = 80, trace = TRUE, B_true = sim_diag$true$B)
fit_short$stop_reason
fit_short$iterations
fit_long <- fit_pfa_wbonly_warm_start(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2,
                                      max_iter = 500, trace = TRUE, B_true = sim_diag$true$B,
                                      tol_mu = 1e-8, tol_phi = 1e-8, also_require_mu_phi = TRUE)
fit_long$stop_reason
fit_long$iterations
# compare_fitted_pi(fit_short, fit_long, sim_diag$X, sim_diag$group)

# std error

eta_ref <- sweep(sim_diag$X %*% fit_long$phi, 2, fit_long$mu, "+") + t(fit_long$lambda_hat[, sim_diag$group, drop = FALSE])
pi_ref <- row_softmax(eta_ref)
sampling_se <- sqrt(pi_ref * (1 - pi_ref) / sim_diag$M)

leak <- compare_fitted_pi(fit_short, fit_long, sim_diag$X, sim_diag$group)
leak
# summary(sampling_se / leak)

leak_ratio <- sampling_se / leak
min(leak_ratio)
which(leak_ratio == min(leak_ratio), arr.ind = TRUE)
mean(leak_ratio)
median(leak_ratio)

# Try Q = 200, face serious problem: leak > noise
sum(leak_ratio < 1)

tail(fit_long$mu_shift_trace, 5)
tail(fit_long$mu_remaining_trace, 5)
tail(fit_long$phi_remaining_trace, 5)
fit_long$tol_mu
fit_long$tol_phi

fit_long <- fit_pfa_wbonly_warm_start(sim_diag$Y, sim_diag$X, sim_diag$group, K = 2,
                                      max_iter = 2000, trace = TRUE, B_true = sim_diag$true$B,
                                      tol = 1e-12, stall_patience = 1e6)
fit_long$stop_reason
fit_long$iterations
tail(fit_long$mu_remaining_trace, 5)

fit_long$is_monotone
fit_long$first_decrease_iter
fit_long$best_iter
tail(diff(fit_long$log_evidence), 20)

tail(fit_long$B_shift_trace, 30)
tail(fit_long$sigma2_shift_trace, 30)
tail(fit_long$B_dist_trace, 30)
which.min(fit_long$B_dist_trace)

# Bdist isn't monotone, check the decreasing rate
fit_long$tol_B
plot(fit_long$B_dist_trace[1:300], type = "l")
abline(v = 28, col = "blue", lty = 2)

first_stable <- which(fit_long$B_shift_trace[1:300] < fit_long$tol_B)[1]
first_stable
fit_long$B_dist_trace[first_stable]

  # Compare wb itsef: warm_start and without warm_start, compare the runtime, accuracy and converged issue
  res_grid <- run_grid(Js = c(50),
                       Qs = c(50),
                       seeds = 2:5,
                       K = 2, 
                       Nj = 15, 
                       sigma2_true = 0.3,
                       M_range = c(50, 200),
                       max_iter = 100, 
                       tol = 1e-4,
                       tol_B = 1e-3,
                       tol_sigma2 = 1e-3,
                       methods = c("glmmTMB")
                       )

dat <- simulate_pfa_data(Q = 30, K = 2, J = 25, N_per_group_range = c(15, 15), P = 3,
                         sigma2 = 0.3, M_range = c(15, 50), seed = 1)

f_cold <- fit_pfa_wbonly(dat$Y, dat$X, dat$group, K = 2, M = dat$M,
                         max_iter = 100, tol = 1e-4, sigma2_init = 0.3,
                         trace = TRUE, B_true = dat$true$B)
f_warm <- fit_pfa_wbonly_warm_start(dat$Y, dat$X, dat$group, K = 2, M = dat$M,
                                    max_iter = 100, tol = 1e-4, sigma2_init = 0.3,
                                    trace = TRUE, B_true = dat$true$B)

f_cold$iterations
f_warm$iterations
mean(f_cold$inner_iter_mean_trace)
mean(f_warm$inner_iter_mean_trace)
plot(f_cold$inner_iter_mean_trace, type = "l", ylim = range(c(f_cold$inner_iter_mean_trace, f_warm$inner_iter_mean_trace)))
lines(f_warm$inner_iter_mean_trace, col = "red")

dat <- simulate_pfa_data(Q = 30, K = 2, J = 25, N_per_group_range = c(15, 16), P = 3,
                         sigma2 = 0.3, M_range = c(15, 50), seed = 1)
fit <- fit_pfa_wbonly_warm_start(dat$Y, dat$X, dat$group, K = 2, M = dat$M,
                                 max_iter = 80, tol = 1e-4,
                                 sigma2_init = 0.3, verbose = FALSE, also_require_mu_phi = FALSE,
                                 trace = TRUE, B_true = dat$true$B)
saveRDS(fit, file = "test_fit.rds")

res_test <- readRDS("test_fit.rds")
res_test
summary(res_test)
