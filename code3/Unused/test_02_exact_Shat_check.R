# Compare the accuracy and runtime of multinomial and PME in Hessian step (in E-step, after Newton iteration for mode)
# One table and one plot output. 1. Table for quantitatively numerical analysis and checking the converged issue for Laplace-EM
# 2. Plot for the calculation rate and acceleration between multinomial and PME, also add the accuracy as box-plot. Two different y-axis

# Finding: 
# 1. NOT CONVERGED, all achieve the maximal of iteration
# 2. Total runtime does not differ too much, however, the Shat step works
# 3. Estimation of B is better when exact_Shat = FALSE attended, need to check the reason

library(bench) # timing
library(tidyverse)
library(ggplot2)

source("sim_data.R")
source("pfa_base.R")
source("pfa_woodbury.R")
source("basic_functions.R")
source("pfa_woodbury_traced.R")

# Parameters and grid setup
Q <- c(50)
seeds <- 1:2
J <- c(50)
K <- 2
N_per_grp <- 15
P <- 3
sigma2 <- 0.3
M_rate <- 150
max_iter <- 100 # Large number checking the converged issue
ll_tol <- 1e-2 # Small tolerance checking the converged issue
sigma2_tol <- 1e-6

all_perms <- function(K) {
  if (K == 1) return(matrix(1, 1, 1))
  sub <- all_perms(K - 1)
  do.call(rbind, lapply(1:K, function(i) cbind(i, sub + (sub >= i))))
}

tucker_congruence <- function(B_hat, B_true) {
  K <- ncol(B_true)
  perms <- all_perms(K)
  best_score <- -Inf
  best_phi <- rep(NA, K)
  for (r in 1:nrow(perms)) {
    p <- perms[r, ]
    phi_k <- sapply(1:K, function(k) {
      a <- B_hat[, k]; b <- B_true[, p[k]]
      sum(a * b) / sqrt(sum(a^2) * sum(b^2))
    })
    score <- sum(abs(phi_k))
    if (score > best_score) { best_score <- score; best_phi <- phi_k }
  }
  list(phi_per_factor = best_phi, phi_mean = mean(abs(best_phi)))
}

rv_coefficient <- function(B_hat, B_true) {
  A <- tcrossprod(B_hat)
  Bm <- tcrossprod(B_true)
  sum(A * Bm) / sqrt(sum(A * A) * sum(Bm * Bm))
}

# Single E-step
shat_fidelity_check <- function(dat, mu, phi, B, sigma2) {
  es_exact <- estep_wb(dat$J, dat$group, dat$Y, dat$X, dat$M, mu, phi, B, sigma2, exact_Shat = TRUE)
  es_surr <- estep_wb(dat$J, dat$group, dat$Y, dat$X, dat$M, mu, phi, B, sigma2, exact_Shat = FALSE)
  rel_fro <- sapply(1:dat$J, function(j) {
    norm(es_surr$S_hat[[j]] - es_exact$S_hat[[j]], "F") / norm(es_exact$S_hat[[j]], "F")
  })
  logdet_exact <- sapply(1:dat$J, function(j) as.numeric(determinant(es_exact$S_hat[[j]], logarithm = TRUE)$modulus))
  logdet_surr <- sapply(1:dat$J, function(j) as.numeric(determinant(es_surr$S_hat[[j]], logarithm = TRUE)$modulus))
  data.frame(group = 1:dat$J, rel_frobenius = rel_fro,
             logdet_exact = logdet_exact, logdet_surr = logdet_surr,
             logdet_abs_diff = abs(logdet_exact - logdet_surr))
}

# Full iterations (very small tolerance) to diagnosis the converged issue
run_paired_arm <- function(dat, exact_Shat, max_iter, tol) {
  t0 <- proc.time()[3]
  fit <- fit_pfa_woodbury_traced(dat$Y, dat$X, dat$group, K = dat$K, M = dat$M,
                                 max_iter = max_iter, tol = tol,
                                 verbose = FALSE, exact_Shat = exact_Shat,
                                 B_true = dat$true$B)
  runtime <- proc.time()[3] - t0
  tc <- tucker_congruence(fit$B, dat$true$B)
  rv <- rv_coefficient(fit$B, dat$true$B)
  summary_row <- data.frame(
    exact_Shat = exact_Shat,
    converged = fit$converged,
    iterations = fit$iterations,
    runtime_sec = as.numeric(runtime),
    log_ev_final = tail(fit$log_evidence, 1),
    sigma2_final = fit$sigma2,
    sigma2_true = dat$true$sigma2,
    sigma2_rel_err = abs(fit$sigma2 - dat$true$sigma2) / dat$true$sigma2,
    B_dist_final = tail(fit$B_dist_trace, 1),
    B_dist_min = min(fit$B_dist_trace),
    tucker_mean = tc$phi_mean,
    rv_coefficient = rv
  )
  list(summary = summary_row, fit = fit)
}

converged_ok <- function(fit, ll_tol, sigma2_tol) {
  le <- fit$log_evidence
  s2 <- fit$sigma2_trace
  n <- length(le)
  if (n < 5) return(FALSE)
  ll_ok <- abs(le[n] - le[n - 1]) < ll_tol
  if (is.null(s2)) return(ll_ok)
  m <- length(s2)
  s2_ok <- abs(s2[m] - s2[m - 1]) < sigma2_tol
  ll_ok && s2_ok
}

first_convergence_iter <- function(fit, ll_tol, sigma2_tol, min_iter = 5) {
  le <- fit$log_evidence
  s2 <- fit$sigma2_trace
  n <- length(le)
  if (n < min_iter)
    return(NA)
  for (i in min_iter:n) {
    ll_ok <- abs(le[i] - le[i - 1]) < ll_tol
    s2_ok <- if (is.null(s2)) TRUE else abs(s2[i] - s2[i - 1]) < sigma2_tol
    if (ll_ok && s2_ok)
      return(i)
  }
  NA
}

setup_shat_inputs <- function(Q, K, N_j, sigma2, seed) {
  set.seed(seed)
  B <- apply_PLT(matrix(rnorm(Q * K, 0, 0.6), Q, K))
  BtB <- crossprod(B)
  M_K <- diag(K) + BtB / sigma2
  M_K_inv <- solve(M_K)
  log_det_MK <- 2 * sum(log(diag(chol(M_K))))
  eta <- matrix(rnorm(N_j * Q), N_j, Q)
  pi <- row_softmax(eta)
  M_j <- rpois(N_j, 150) + 1
  Mpi <- sweep(pi, 1, M_j, "*")
  d <- as.numeric(colSums(Mpi))
  list(B = B, M_K = M_K, M_K_inv = M_K_inv, log_det_MK = log_det_MK,
       pi = pi, Mpi = Mpi, d = d, sigma2 = sigma2, Q = Q, K = K)
}

shat_exact_step <- function(inp) {
  Sinv_mat <- diag(1 / inp$sigma2, inp$Q) - inp$B %*% inp$M_K_inv %*% t(inp$B) / inp$sigma2^2
  negH <- diag(inp$d, inp$Q) - crossprod(inp$pi, inp$Mpi) + Sinv_mat
  cL  <- chol(negH)
  S_hat <- chol2inv(cL)
  log_det_shat <- -2 * sum(log(diag(cL)))
  invisible(NULL)
}

shat_surrogate_step <- function(inp) {
  dt  <- inp$d + 1 / inp$sigma2
  idt <- 1 / dt
  Gam <- t(inp$B) %*% (idt * inp$B) - inp$sigma2^2 * inp$M_K
  L <- chol(-Gam + 1e-10 * diag(inp$K))
  log_det_shat <- -(sum(log(dt)) - 2 * inp$K * log(inp$sigma2) - inp$log_det_MK +
                      2 * sum(log(diag(L))))
  GamInv <- chol2inv(L)
  IdtB <- idt * inp$B
  S_hat  <- diag(idt, inp$Q) + IdtB %*% GamInv %*% t(IdtB)
  invisible(NULL)
}

benchmark_shat_scaling <- function(Q_grid, K, N_j, sigma2, seed) {
  results <- lapply(Q_grid, function(Q) {
    inp <- setup_shat_inputs(Q, K, N_j, sigma2, seed + Q)
    mk <- bench::mark(
      exact = shat_exact_step(inp),
      surrogate = shat_surrogate_step(inp),
      check = FALSE, min_iterations = 5
    )
    data.frame(
      Q = Q, K = K,
      time_exact_ms = as.numeric(mk$median[1]) * 1000,
      time_surrogate_ms = as.numeric(mk$median[2]) * 1000,
      speedup_observed = as.numeric(mk$median[1]) / as.numeric(mk$median[2]),
      speedup_theoretical_Q_over_K = Q / K
    )
  })
  do.call(rbind, results)
}

# Plot + table: pure Hessian-step timing, computed once up front so every (q, j, seed) cell below can look up its own row by q
Q_grid <- sort(unique(c(25, 50, 100, 200, 400, 800, 1600, Q)))
shat_timing <- benchmark_shat_scaling(Q_grid, K, N_per_grp, sigma2, seeds[1])
rownames(shat_timing) <- NULL

shat_plot <- shat_timing %>%
  pivot_longer(cols = c(time_exact_ms, time_surrogate_ms),
               names_to = "branch", values_to = "time_ms") %>%
  mutate(branch = recode(branch,
                         time_exact_ms = "exact_Shat = TRUE",
                         time_surrogate_ms = "exact_Shat = FALSE")) %>%
  ggplot(aes(x = Q, y = time_ms, color = branch)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "S_hat construction time vs Q", x = "Q", y = "time per call (ms)", color = NULL) +
  theme_minimal()

print(shat_plot)
print(shat_timing)

# Run the code on grid
all_results <- list()
all_snapshots <- list()

for (q in Q) {
  for (j in J) {
    for (seed in seeds) {
      
      dat <- simulate_pfa_data(Q = q, K = K, J = j, N_per_group = N_per_grp,
                               P = P, sigma2 = sigma2, M_rate = M_rate, seed = seed)
      
      # Use true value of theta
      snapshot <- shat_fidelity_check(dat, dat$true$mu, dat$true$phi, dat$true$B, dat$true$sigma2)
      
      arm_exact <- run_paired_arm(dat, exact_Shat = TRUE, max_iter, tol = 1e-10)
      arm_surr  <- run_paired_arm(dat, exact_Shat = FALSE, max_iter, tol = 1e-10)
      
      results <- rbind(arm_exact$summary, arm_surr$summary)
      results$converged_at_max_iter <- c(
        converged_ok(arm_exact$fit, ll_tol, sigma2_tol),
        converged_ok(arm_surr$fit, ll_tol, sigma2_tol)
      )
      results$first_convergence_iter <- c(
        first_convergence_iter(arm_exact$fit, ll_tol, sigma2_tol),
        first_convergence_iter(arm_surr$fit, ll_tol, sigma2_tol)
      )
      
      shat_row <- shat_timing[shat_timing$Q == q, ]
      results$shat_step_ms <- c(shat_row$time_exact_ms, shat_row$time_surrogate_ms)
      
      results$Q <- q
      results$J <- j
      results$seed <- seed
      
      cat(sprintf("Q=%d J=%d seed=%d\n", q, j, seed))
      print(summary(snapshot$rel_frobenius))
      print(summary(snapshot$logdet_abs_diff))
      print(results)
      
      all_results[[length(all_results) + 1]] <- results
      all_snapshots[[length(all_snapshots) + 1]] <- cbind(Q = q, J = j, seed = seed, snapshot)
    }
  }
}

results_all <- do.call(rbind, all_results)
snapshot_all <- do.call(rbind, all_snapshots)
rownames(results_all) <- NULL
rownames(snapshot_all) <- NULL

print(results_all)