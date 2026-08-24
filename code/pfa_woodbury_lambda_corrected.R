correct_lambda_edgeworth <- function(lambda_hat, S_hat, Y, X, group, M, mu, phi) {
  Q <- nrow(lambda_hat); J <- ncol(lambda_hat)
  stopifnot(length(S_hat) == J)

  lam_corr <- lambda_hat
  rel_size <- numeric(J)
  Xphi <- X %*% phi

  for (j in seq_len(J)) {
    idx <- which(group == j)
    n_j <- length(idx)
    Sj  <- S_hat[[j]]
    diagS <- diag(Sj)

    eta <- matrix(mu, n_j, Q, byrow = TRUE) +
           Xphi[idx, , drop = FALSE] +
           matrix(lambda_hat[, j], n_j, Q, byrow = TRUE)
    eta <- eta - apply(eta, 1, max)   
    pi_mat <- exp(eta); pi_mat <- pi_mat / rowSums(pi_mat)  

    TS <- numeric(Q)
    for (i in seq_len(n_j)) {
      pi_i  <- pi_mat[i, ]
      Spi_i <- as.vector(Sj %*% pi_i)
      piDiagS <- sum(pi_i * diagS)
      piSpi   <- sum(pi_i * Spi_i)
      TS <- TS + M[idx[i]] * pi_i * (diagS - 2 * Spi_i - piDiagS + 2 * piSpi)
    }

    mu1 <- -0.5 * as.vector(Sj %*% TS)

    lam_corr[, j] <- lambda_hat[, j] + mu1
    rel_size[j] <- sqrt(sum(mu1^2)) / max(sqrt(sum(lambda_hat[, j]^2)), 1e-8)
  }

  list(lambda_corrected = lam_corr, rel_size = rel_size)
}

g_map_corrected <- function(sigma2_prior, dat, J_fixed, K_true, correct = FALSE,
                             estep_max_iter = 200, estep_gtol = 1e-6) {
  Q <- ncol(dat$Y)
  mu  <- rep(0, Q)
  phi <- matrix(0, ncol(dat$X), Q)

  es <- estep_wb(J_fixed, dat$group, dat$Y, dat$X, dat$M, mu, phi,
                 dat$true$B, sigma2_prior,
                 max_iter = estep_max_iter, gtol = estep_gtol)

  if (is.null(es$lambda_hat)) {
    stop("estep_wb() output has no $lambda_hat - run names(es) and check what changed.")
  }
  lam <- es$lambda_hat

  rel_size <- rep(NA_real_, J_fixed)
  if (correct) {
    cc <- correct_lambda_edgeworth(lam, es$S_hat, dat$Y, dat$X, dat$group, dat$M, mu, phi)
    lam <- cc$lambda_corrected
    rel_size <- cc$rel_size
  }

  S_obs <- tcrossprod(lam) / J_fixed
  for (j in 1:J_fixed) S_obs <- S_obs + es$S_hat[[j]] / J_fixed

  eigs <- eigen(S_obs, symmetric = TRUE, only.values = TRUE)$values
  sigma2_out <- mean(sort(eigs)[1:(Q - K_true)])
  
  rt <- rubin_thayer_wb(S_obs, K_true, B_init = dat$true$B, sigma2_init = sigma2_prior)
  B_out  <- apply_PLT(rt$B)
  B_dist <- subspace_dist(B_out, dat$true$B)

  list(sigma2_out = sigma2_out, rel_size = rel_size, B_out = B_out, B_dist = B_dist)
}

run_single_round_comparison <- function(use_test16_seed = TRUE) {

  N_per_grp <- 15
  K_true <- 2
  J_fixed <- 6
  Q_grid <- c(50)
  target <- (J_fixed - K_true) / J_fixed

  cat(sprintf("%5s | %10s | %10s %10s | %10s %10s | %10s\n",
              "Q", "s2_true", "s2_base", "s2_corr", "Bd_base", "Bd_corr", "corr_sz"))
  for (Qv in Q_grid) {
    seed <- if (use_test16_seed) 2000 * Qv + 1 else sample.int(1e6, 1)
    
    dat <- simulate_pfa_data(Q = Qv, K = K_true, J = J_fixed,
                             N_per_group = N_per_grp, seed = seed, sigma2 = target)
    if (abs(dat$true$sigma2 - target) > 1e-3) {
      stop(sprintf("dat$true$sigma2 (%.4f) != target (%.4f) - fix before trusting anything below.",
                   dat$true$sigma2, target))
    }
    cat(sprintf("  dat$true$sigma2 = %.4f matches target = %.4f (OK)\n", dat$true$sigma2, target))

    base <- g_map_corrected(target, dat, J_fixed, K_true, correct = FALSE)
    corr <- g_map_corrected(target, dat, J_fixed, K_true, correct = TRUE)

    cat(sprintf("%5d | %10.4f | %10.4f %10.4f | %10.4f %10.4f | %10.4f\n",
                Qv, dat$true$sigma2, base$sigma2_out, corr$sigma2_out,
                base$B_dist, corr$B_dist,
                mean(corr$rel_size, na.rm = TRUE)))
  }
}