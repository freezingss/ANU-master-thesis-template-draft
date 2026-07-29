# =============================================================================
# test_02_woodbury_validation.R
# =============================================================================

source("pfa_base_dense_armijo.R")
source("pfa_woodbury_armijo.R") 
source("sim_data.R")

set.seed(3)
tolm <- 1e-8

Q <- 40; K <- 2; Nj <- 20
dat <- simulate_pfa_data(Q = Q, K = K, J = 1, N_per_group = Nj, seed = 3)
Y_j <- dat$Y; X_j <- dat$X; M_j <- dat$M
tr <- dat$true
B <- tr$B; s2 <- tr$sigma2; mu <- tr$mu; phi <- tr$phi

BtB <- crossprod(B)
M_K <- diag(K) + BtB / s2
M_K_inv <- solve(M_K)
ldMK <- 2 * sum(log(diag(chol(M_K))))
Sigma <- B %*% t(B) + s2 * diag(Q)

check <- function(name, err, tol = tolm) {
  cat(sprintf("  %-52s %s   (err = %.2e)\n", name,
              ifelse(err < tol, "PASS", "FAIL"), err))
}

cat("=============================================================\n")
cat(" (a)-(c) Woodbury linear algebra exactness\n")
cat("=============================================================\n")

# (a) Sigma^{-1} v
v <- rnorm(Q)
wb_Sv  <- v / s2 - B %*% (M_K_inv %*% (t(B) %*% v)) / s2^2
check("(a) Sigma^{-1} v : Woodbury vs dense solve",
      max(abs(wb_Sv - solve(Sigma, v))))

# (b) DLR surrogate Newton direction
fixed <- matrix(rep(mu, each = Nj), Nj, Q) + X_j %*% phi
lambda <- rnorm(Q, 0, 0.3)
eta <- sweep(fixed, 2, lambda, "+")
pi <- row_softmax(eta); Mpi <- sweep(pi, 1, M_j, "*")
d <- as.numeric(colSums(Mpi))
g <- rnorm(Q)

dt <- d + 1 / s2; idt <- 1 / dt
A <- diag(dt, Q) - B %*% M_K_inv %*% t(B) / s2^2      # DLR surrogate -Hess
Gam <- t(B) %*% (idt * B) - s2^2 * M_K
L <- chol(-Gam)
wb_dir <- idt * g - idt * (B %*% backsolve(L, forwardsolve(t(L),
                                                           -crossprod(B, idt * g))))
check("(b) A^{-1} g : Woodbury vs dense solve", max(abs(wb_dir - solve(A, g))))

# (c) Level-2 beta identity
check("(c) Sigma^{-1} B == (1/sigma2) B M_K^{-1}",
      max(abs(solve(Sigma, B) - B %*% M_K_inv / s2)))

# (d) Level-3 analytic log-determinant
cat("=============================================================\n")
cat(" (d) Level-3 analytic log|S_hat|\n")
cat("=============================================================\n")
ld_analytic <- -(sum(log(dt)) - 2 * K * log(s2) - ldMK + 2 * sum(log(diag(L))))
ld_dense <- -as.numeric(determinant(A, logarithm = TRUE)$modulus)
check("(d) log|S_hat| : analytic (O(K^3)) vs dense (O(Q^3))",
      abs(ld_analytic - ld_dense))

cat("=============================================================\n")
cat(" (e) Stopping-rule defect and identical modes\n")
cat("=============================================================\n")

Sigma_inv <- solve(Sigma)
r_dense <- laplace_lambda_j_dense(Y_j, X_j, M_j, mu, phi, Sigma_inv,
                                 max_iter = 200, gtol = 1e-8)

lap_wb_stepstop <- function(max_iter = 50, tol = 1e-9) {
  Sinv <- function(v) v / s2 - B %*% (M_K_inv %*% (t(B) %*% v)) / s2^2
  lp <- function(a) { e <- sweep(fixed, 2, a, "+")
    sum(Y_j * e) - sum(M_j * row_logsumexp(e)) - 0.5 * sum(a * Sinv(a)) }
  a <- rep(0, Q); lc <- lp(a); gmax <- Inf
  for (it in 1:max_iter) {
    e <- sweep(fixed, 2, a, "+"); p <- row_softmax(e)
    Mp <- sweep(p, 1, M_j, "*")
    dv <- as.numeric(colSums(Mp)); gg <- as.numeric(colSums(Y_j - Mp)) - Sinv(a)
    gmax <- max(abs(gg))
    dtv <- dv + 1 / s2; idtv <- 1 / dtv
    Gm <- t(B) %*% (idtv * B) - s2^2 * M_K
    Lc <- chol(-Gm)
    dir <- idtv * gg - idtv * (B %*% backsolve(Lc, forwardsolve(t(Lc),
                                                                -crossprod(B, idtv * gg))))
    st <- 1
    for (bt in 1:20) { if (lp(a + st * dir) >= lc - 1e-12) break; st <- st / 2 }
    a <- a + st * dir; lc <- lp(a)
    if (max(abs(st * dir)) < tol) break        # <- the defective test
  }
  list(lambda = a, iters = it, grad = gmax)
}
r_old <- lap_wb_stepstop()
cat(sprintf("  step-size stopping: 'converged' at iter %d with max|g| = %.3f",
            r_old$iters, r_old$grad))
cat("  <- FALSE convergence\n")

r_wb <- laplace_lambda_j_wb2(Y_j, X_j, M_j, mu, phi, B, s2,
                            BtB, M_K, M_K_inv, ldMK,
                            max_iter = 400, gtol = 1e-6)
cat(sprintf("  gradient stopping : iter %d, max|g| = %.2e\n",
            r_wb$n_iter, r_wb$grad_norm))
cat(sprintf("  cor(lambda_dense, lambda_wb) = %.6f   max|diff| = %.2e\n",
            cor(r_dense$lambda_hat, r_wb$lambda_hat),
            max(abs(r_dense$lambda_hat - r_wb$lambda_hat))))
check("(e) identical posterior modes (dense vs Woodbury)",
      max(abs(r_dense$lambda_hat - r_wb$lambda_hat)), tol = 1e-3)

cat("=============================================================\n")
cat(" (f) Per-Newton-step microbenchmark (this machine)\n")
cat("=============================================================\n")
cat(sprintf("%6s %12s %12s %10s %12s\n",
            "Q", "dense(ms)", "wb(ms)", "speedup", "lp()(ms)"))
for (Qb in c(50, 100, 200, 400)) {
  db <- simulate_pfa_data(Q = Qb, K = 2, J = 1, N_per_group = 20, seed = Qb)
  trb <- db$true; Bb <- trb$B; sb <- trb$sigma2
  MKb <- diag(2) + crossprod(Bb) / sb
  MKbi <- solve(MKb)
  Sig <- Bb %*% t(Bb) + sb * diag(Qb)
  Sigi <- solve(Sig)
  fx <- matrix(rep(trb$mu, each = 20), 20, Qb) + db$X %*% trb$phi
  pb <- row_softmax(fx); Mpb <- sweep(pb, 1, db$M, "*")
  dv <- as.numeric(colSums(Mpb)); gb <- rnorm(Qb)
  Hb <- diag(dv, Qb) - crossprod(pb, Mpb) + Sigi
  nb <- 500                        # many reps: the Woodbury solve is sub-ms
  t_d <- system.time(for (r in 1:nb) x1 <- solve(Hb, gb))["elapsed"] / nb * 1e3
  dtb <- dv + 1 / sb; idtb <- 1 / dtb
  t_w <- system.time(for (r in 1:nb) {
    Gm <- t(Bb) %*% (idtb * Bb) - sb^2 * MKb
    Lb <- chol(-Gm)
    x2 <- idtb * gb - idtb * (Bb %*% backsolve(Lb, forwardsolve(t(Lb),
                                                                -crossprod(Bb, idtb * gb))))
  })["elapsed"] / nb * 1e3
  lpf <- function(a) { e <- sweep(fx, 2, a, "+")
    sum(db$Y * e) - sum(db$M * row_logsumexp(e)) }
  a0 <- rnorm(Qb, 0, 0.1)
  t_l <- system.time(for (r in 1:nb) z <- lpf(a0))["elapsed"] / nb * 1e3
  t_w <- max(t_w, 1e-3)            # floor at clock resolution for the ratio
  cat(sprintf("%6d %12.3f %12.3f %9.1fx %12.3f\n", Qb, t_d, t_w, t_d / t_w, t_l))
}
cat("Note: the lp() column is the cost of ONE backtracking evaluation.  A\n")
cat("Newton step performs ~1-5 of them in BOTH methods, so in interpreted R\n")
cat("lp() becomes the dominant per-step cost once the Woodbury solve has\n")
cat("collapsed to ~0.05 ms; the full O(Q^3 -> QK^2) gain appears when lp()\n")
cat("is also compiled (C++/Julia).  Full-scale numbers in the file header.\n\n")

cat("=============================================================\n")
cat(" (g) Full-EM smoke test: dense vs Woodbury (surrogate & exact S_hat)\n")
cat("=============================================================\n")
dat2 <- simulate_pfa_data(Q = 30, K = 2, J = 20, N_per_group = 12, seed = 7)

t_fd <- system.time(
  f_d <- fit_pfa_dense(dat2$Y, dat2$X, dat2$group, K = 2, M = dat2$M,
                       max_iter = 15, verbose = FALSE))["elapsed"]

# Woodbury with the surrogate (Poisson) S_hat -- production / fast path
t_fw_s <- system.time(
  f_w_s <- fit_pfa_woodbury(dat2$Y, dat2$X, dat2$group, K = 2, M = dat2$M,
                            max_iter = 15, verbose = FALSE,
                            exact_Shat = FALSE))["elapsed"]

# Woodbury with the exact multinomial S_hat -- should match dense
t_fw_e <- system.time(
  f_w_e <- fit_pfa_woodbury(dat2$Y, dat2$X, dat2$group, K = 2, M = dat2$M,
                            max_iter = 15, verbose = FALSE,
                            exact_Shat = TRUE))["elapsed"]

cat(sprintf("  dense       : %2d iters, %5.1fs, final log_ev = %.3f\n",
            f_d$iterations,   t_fd,   tail(f_d$log_evidence, 1)))
cat(sprintf("  wb-surrogate : %2d iters, %5.1fs, final log_ev = %.3f\n",
            f_w_s$iterations, t_fw_s, tail(f_w_s$log_evidence, 1)))
cat(sprintf("  wb-exact     : %2d iters, %5.1fs, final log_ev = %.3f\n",
            f_w_e$iterations, t_fw_e, tail(f_w_e$log_evidence, 1)))

cat(sprintf("  sigma2   : dense %.4f | wb-surr %.4f | wb-exact %.4f | true %.3f\n",
            f_d$sigma2, f_w_s$sigma2, f_w_e$sigma2, dat2$true$sigma2))
cat(sprintf("  B subspace dist : dense %.4f | wb-surr %.4f | wb-exact %.4f\n",
            subspace_dist(f_d$B,   dat2$true$B),
            subspace_dist(f_w_s$B, dat2$true$B),
            subspace_dist(f_w_e$B, dat2$true$B)))
cat(sprintf("  lambda mean_r   : dense %.4f | wb-surr %.4f | wb-exact %.4f\n",
            lambda_mean_cor(f_d$lambda_hat,   dat2$true$lambda),
            lambda_mean_cor(f_w_s$lambda_hat, dat2$true$lambda),
            lambda_mean_cor(f_w_e$lambda_hat, dat2$true$lambda)))

mono_d <- all(diff(f_d$log_evidence)   > -1e-4 * abs(f_d$log_evidence[-1]))
mono_w_s <- all(diff(f_w_s$log_evidence) > -1e-4 * abs(f_w_s$log_evidence[-1]))
mono_w_e <- all(diff(f_w_e$log_evidence) > -1e-4 * abs(f_w_e$log_evidence[-1]))
cat(sprintf("  log-evidence non-decreasing: dense %s | wb-surr %s | wb-exact %s\n",
            mono_d, mono_w_s, mono_w_e))

cat(sprintf("  |log_ev  dense - wb-exact| = %.4f   (~0 expected: same true S_hat)\n",
            abs(tail(f_d$log_evidence, 1) - tail(f_w_e$log_evidence, 1))))
cat(sprintf("  |log_ev  dense - wb-surr| = %.4f    (nonzero: surrogate S_hat bias)\n",
            abs(tail(f_d$log_evidence, 1) - tail(f_w_s$log_evidence, 1))))
cat(sprintf("  |sigma2  dense - wb-exact| = %.4f   (~0 expected)\n",
            abs(f_d$sigma2 - f_w_e$sigma2)))
cat("  (wb-surrogate uses the Poisson S_hat, which under-states posterior\n")
cat("   variance -- Prop 5.4.1; wb-exact uses the true multinomial S_hat and\n")
cat("   should track dense to numerical tolerance.)\n")