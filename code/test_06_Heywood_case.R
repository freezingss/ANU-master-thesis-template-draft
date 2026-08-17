# test_06_Heywood_case
# -> Check the classic target (J-K)/J
# -> Whether exact_Shat influences the estimation?
# -> Does growing J alongside Q rescue sigma2 recovery?

source("pfa_woodbury.R")
source("sim_data.R")

g_map <- function(sigma2_prior, dat, J_fixed, K_true, exact_Shat = FALSE) {
  Q <- ncol(dat$Y)
  mu <- rep(0, Q); phi <- matrix(0, ncol(dat$X), Q)
  es <- estep_wb(J_fixed, dat$group, dat$Y, dat$X, dat$M, mu, phi,
                 dat$true$B, sigma2_prior, max_iter = 200, gtol = 1e-6,
                 exact_Shat = exact_Shat)
  S_obs <- tcrossprod(es$lambda_hat) / J_fixed
  for (j in 1:J_fixed) S_obs <- S_obs + es$S_hat[[j]] / J_fixed
  eigs <- eigen(S_obs, symmetric = TRUE, only.values = TRUE)$values
  mean(sort(eigs)[1:(Q - K_true)])
}

find_fixed_point <- function(dat, J_fixed, K_true, lo = 0.02, hi = 2.0, tol = 1e-3,
                             exact_Shat = FALSE) {
  f <- function(x) g_map(x, dat, J_fixed, K_true, exact_Shat = exact_Shat) - x
  flo <- f(lo); fhi <- f(hi)
  if (sign(flo) == sign(fhi)) return(NA_real_)
  for (it in 1:40) {
    mid <- (lo + hi) / 2
    fmid <- f(mid)
    if (abs(fmid) < tol || (hi - lo) < tol) return(mid)
    if (sign(fmid) == sign(flo)) { lo <- mid; flo <- fmid } else { hi <- mid }
  }
  (lo + hi) / 2
}

J_fixed <- 6
K_true <- 2
target <- (J_fixed - K_true) / J_fixed

part1_results <- list()

for (Q in c(50, 400)) {
  cat(sprintf("Part 1: J=%d fixed, Q = %d \n", J_fixed, Q))
  dat <- simulate_pfa_data(Q = Q, K = K_true, J = J_fixed, N_per_group = 15,
                           sigma2 = target, seed = 2000 * Q + 1)
  grid <- c(0.02, 0.05, 0.1, 0.2, 0.3, target, 0.5, 0.75, 1.0, 1.5, 2.0)
  cat(sprintf("%12s | %12s | %10s\n", "sigma2_prior", "g(sigma2_prior)", "g > input?"))
  cat(strrep("-", 45), "\n")
  gvals <- numeric(length(grid))
  for (i in seq_along(grid)) {
    gvals[i] <- g_map(grid[i], dat, J_fixed, K_true, exact_Shat = FALSE)
    cat(sprintf("%12.4f | %12.4f | %10s\n",
                grid[i], gvals[i], ifelse(gvals[i] > grid[i], "YES (grows)", "no (shrinks)")))
  }
  diffs <- gvals - grid
  sign_change <- which(diff(sign(diffs)) != 0)
  if (length(sign_change) > 0) {
    cat(sprintf("\n  Sign change (approx fixed point) between sigma2_prior = %.3f and %.3f\n",
                grid[sign_change[1]], grid[sign_change[1] + 1]))
  } else {
    cat("\n  No sign change detected in g(x)-x over this grid.\n")
  }
  g_at_target <- gvals[grid == target]   # already computed above, no need to re-call g_map
  cat(sprintf("\n  Reference: g(true sigma2=%.4f) = %.4f  (target %.4f)\n\n",
              target, g_at_target, target))
  
  part1_results[[paste0("Q", Q)]] <- data.frame(sigma2_prior = grid, g = gvals)
  saveRDS(part1_results, "test_06_part1_results.rds")
}

K_true <- 2
Q_fixed <- 400
J_grid <- c(6, 10, 20, 40, 80)

cat(sprintf("Part 2: Q=%d fixed, J-sweep \n", Q_fixed))
cat(sprintf("%6s | %10s | %14s | %16s\n", "J", "target", "sigma2_fp", "sigma2_fp/target"))
cat(strrep("-", 55), "\n")

part2_results <- data.frame()
for (J_fixed in J_grid) {
  target <- (J_fixed - K_true) / J_fixed
  dat <- simulate_pfa_data(Q = Q_fixed, K = K_true, J = J_fixed, N_per_group = 15,
                           sigma2 = target, seed = 5000 * J_fixed + 1)
  fp <- find_fixed_point(dat, J_fixed, K_true, exact_Shat = FALSE)
  cat(sprintf("%6d | %10.4f | %14.4f | %16.4f\n",
              J_fixed, target, fp, fp / target))
  part2_results <- rbind(part2_results,
                         data.frame(J = J_fixed, target = target,
                                    sigma2_fp = fp, ratio = fp / target))
  saveRDS(part2_results, "test_06_part2_results.rds")
}

K_true <- 2
Q_grid <- c(50, 100, 200, 400, 800)
J_grid3 <- c(10, 20, 40, 80, 160)   

part3_results <- data.frame()

for (Q in Q_grid) {
  cat(sprintf("\nPart 3: Q = %d, J-sweep\n", Q))
  cat(sprintf("%6s | %10s | %14s | %16s\n", "J", "target", "sigma2_fp", "ratio"))
  cat(strrep("-", 55), "\n")
  
  for (J in J_grid3) {
    target <- (J - K_true) / J
    dat <- simulate_pfa_data(Q = Q, K = K_true, J = J, N_per_group = 15,
                             sigma2 = target, seed = 9000 * Q + 17 * J + 1)
    fp <- find_fixed_point(dat, J, K_true, exact_Shat = FALSE)
    ratio <- fp / target
    cat(sprintf("%6d | %10.4f | %14.4f | %16.4f\n", J, target, fp, ratio))
    part3_results <- rbind(part3_results,
                           data.frame(Q = Q, J = J, target = target,
                                      sigma2_fp = fp, ratio = ratio))
  }
  saveRDS(part3_results, "test_06_part3_results.rds")
}



j_needed_for_ratio <- function(df_q, target_ratio = 0.9) {
  df_q <- df_q[order(df_q$J), ]
  df_q <- df_q[!is.na(df_q$ratio), ]
  if (nrow(df_q) < 2) return(NA_real_)
  if (all(df_q$ratio < target_ratio)) return(NA_real_)  
  above <- which(df_q$ratio >= target_ratio)[1]
  if (above == 1) return(df_q$J[1])  
  below <- above - 1
  x0 <- log(df_q$J[below]); x1 <- log(df_q$J[above])
  y0 <- df_q$ratio[below];  y1 <- df_q$ratio[above]
  exp(x0 + (target_ratio - y0) * (x1 - x0) / (y1 - y0))
}

j_needed_summary <- do.call(rbind, lapply(split(part3_results, part3_results$Q),
                                          function(df_q) data.frame(Q = df_q$Q[1], J_needed = j_needed_for_ratio(df_q, 0.9))))
j_needed_summary <- j_needed_summary[order(j_needed_summary$Q), ]
saveRDS(j_needed_summary, "test_06_part3_j_needed.rds")

cat("\nJ needed to reach ratio >= 0.9, by Q:\n")
print(j_needed_summary)

valid <- !is.na(j_needed_summary$J_needed)
if (sum(valid) >= 2) {
  fit <- lm(log(J_needed) ~ log(Q), data = j_needed_summary[valid, ])
  cat(sprintf("\nEmpirical exponent: J_needed ~ Q^%.2f  (theory conjectures ~Q^2, up to log Q)\n",
              coef(fit)[2]))
  saveRDS(fit, "test_06_part3_exponent_fit.rds")
}