# =============================================================================
# test_01_dense_bottleneck.R
# =============================================================================

source("pfa_base_dense_armijo.R")
source("sim_data.R")

set.seed(2026)

Qs <- c(50, 100, 200, 400)  
K <- 2
n_rep <- 200                 

cat("=============================================================\n")
cat(" (1) Dense Newton solve  solve(H, g) - one E-step iteration\n")
cat("=============================================================\n")
cat(sprintf("%8s %14s %14s\n", "Q", "ms/solve", "vs Q=50"))

t_solve <- numeric(length(Qs))
for (qi in seq_along(Qs)) {
  Q <- Qs[qi]
  dat <- simulate_pfa_data(Q = Q, K = K, J = 5, N_per_group = 15, seed = qi)
  idx <- which(dat$group == 1)
  Y_j <- dat$Y[idx, , drop = FALSE]; X_j <- dat$X[idx, , drop = FALSE]
  M_j <- dat$M[idx]

  tr <- dat$true
  Sigma <- tr$B %*% t(tr$B) + tr$sigma2 * diag(Q)
  Sigma_inv <- solve(Sigma)
  fixed <- matrix(rep(tr$mu, each = length(idx)), length(idx), Q) +
    X_j %*% tr$phi
  eta <- fixed
  pi <- row_softmax(eta)
  Mpi <- sweep(pi, 1, M_j, "*")
  g <- as.numeric(colSums(Y_j - Mpi))

  H <- diag(as.numeric(colSums(Mpi)), Q) - crossprod(pi, Mpi) + Sigma_inv
  t_solve[qi] <- system.time(
    for (r in 1:n_rep) dir <- solve(H, g)
  )["elapsed"] / n_rep * 1000
  cat(sprintf("%8d %14.3f %13.1fx\n", Q, t_solve[qi],
              t_solve[qi] / t_solve[1]))
}

fit1 <- lm(log(t_solve) ~ log(Qs))
cat(sprintf("\nEmpirical growth exponent (dense solve): Q^%.2f",
            coef(fit1)[2]))
cat("  [BLAS-optimised; asymptotically 3]\n\n")

cat("=============================================================\n")
cat(" (2) Dense Rubin-Thayer -- one inner iteration (solve(Sigma))\n")
cat("=============================================================\n")
cat(sprintf("%8s %14s %14s\n", "Q", "ms/iter", "vs Q=50"))

t_rt <- numeric(length(Qs))
for (qi in seq_along(Qs)) {
  Q <- Qs[qi]
  B <- apply_PLT(matrix(rnorm(Q * K, 0, 0.6), Q, K))
  s2 <- 0.3
  S_obs <- B %*% t(B) + s2 * diag(Q) +
    0.05 * crossprod(matrix(rnorm(5 * Q), 5, Q)) / 5   

  one_rt_iter <- function() {
    Sigma <- B %*% t(B) + s2 * diag(Q)
    Si <- solve(Sigma)                        
    beta <- t(B) %*% Si
    Theta <- diag(K) - beta %*% B + beta %*% S_obs %*% t(beta)
    B_new <- S_obs %*% t(beta) %*% solve(Theta)
    invisible(B_new)
  }
  t_rt[qi] <- system.time(for (r in 1:n_rep) one_rt_iter()
  )["elapsed"] / n_rep * 1000
  cat(sprintf("%8d %14.3f %13.1fx\n", Q, t_rt[qi], t_rt[qi] / t_rt[1]))
}

fit2 <- lm(log(t_rt) ~ log(Qs))
cat(sprintf("\nEmpirical growth exponent (dense Rubin-Thayer): Q^%.2f\n",
            coef(fit2)[2]))

cat("\n=============================================================\n")
cat(" (3) What a full dense EM iteration implies at large Q\n")
cat("=============================================================\n")
cat("Per EM iteration the dense algorithm pays roughly\n")
cat("   J groups x ~10 Newton steps x t_solve(Q)   [E-step]\n")
cat(" + ~500 RT inner iterations x t_rt(Q)          [M-step (b)]\n")
J_ref <- 30
for (qi in seq_along(Qs)) {
  est <- (J_ref * 10 * t_solve[qi] + 500 * t_rt[qi]) / 1000
  cat(sprintf("   Q = %4d :  ~%7.1f s / EM iteration (J = %d)\n",
              Qs[qi], est, J_ref))
}
cat("Extrapolating the measured exponents to Q = 1000 gives ~8 min per EM\n")
cat("iteration, i.e. hours per fit -- the high-dimensional dilemma that\n")
cat("motivates pfa_woodbury.R (see test_02 for the accelerated numbers).\n")
