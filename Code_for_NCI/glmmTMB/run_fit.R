run_single_fit <- function(J, Q, K, seed, method,
                           N_per_group_range, M_range, P = 3,
                           sigma2_true = 0.3, trace = FALSE, ...) {
  
  dat <- simulate_pfa_data(Q = Q, K = K, J = J, N_per_group_range = N_per_group_range,
                           P = P, sigma2 = sigma2_true, M_range = M_range, seed = seed)
  
  fit <- switch(method,
                wb = fit_pfa_wbonly_warm_start(dat$Y, dat$X, dat$group, K, M = dat$M,
                                               B_true = dat$true$B, sigma2_true = dat$true$sigma2,
                                               trace = trace, ...),
                glmmTMB = fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K),
                stop(sprintf("method '%s' is not yet wired into run_single_fit", method)))
  
  list(config = list(J = J, Q = Q, K = K, seed = seed, method = method,
                     N_per_group_range = N_per_group_range, M_range = M_range,
                     P = P, sigma2_true_input = sigma2_true),
       truth = c(dat$true, list(n_per_group = dat$n_per_group)),
       fit = fit)
}
