# Observation Level Outliers
simulate_pfa_data <- function(Q, K, J, N_per_group_range, P, sigma2, M_range,
                              dispersion = Inf,
                              outlier_frac = 0,
                              outlier_type = c("spike", "random_composition", "shock"),
                              outlier_severity = NULL,
                              seed) {
  
  outlier_type <- match.arg(outlier_type)
  if (!is.null(seed)) set.seed(seed)
  
  mu <- rnorm(Q, 0, 0.4)
  mu <- mu - mean(mu)
  phi <- matrix(0, P, Q)
  if (P > 1) phi[2:P, ] <- rnorm((P - 1) * Q, 0, 0.3)
  B <- apply_PLT(matrix(rnorm(Q * K, 0, 0.6), Q, K))
  Fm <- matrix(rnorm(J * K), J, K)
  lambda <- B %*% t(Fm) + matrix(rnorm(Q * J, 0, sqrt(sigma2)), Q, J)
  
  n_j <- sample(N_per_group_range[1]:N_per_group_range[2], J, replace = TRUE)
  N <- sum(n_j)
  group <- rep(seq_len(J), times = n_j)
  
  X <- cbind(1, matrix(rnorm(N * (P - 1)), N, P - 1))
  M <- sample(M_range[1]:M_range[2], N, replace = TRUE)
  
  eta <- sweep(X %*% phi, 2, mu, "+") + t(lambda[, group, drop = FALSE])
  
  Y <- t(sapply(seq_len(N), function(i) {
    p_raw <- exp(eta[i, ] - max(eta[i, ]))
    p <- p_raw / sum(p_raw)
    if (is.finite(dispersion)) {
      g <- rgamma(Q, shape = p * dispersion, rate = 1)
      p <- g / sum(g)
    }
    as.vector(rmultinom(1, M[i], p))
  }))
  
  n_out <- round(N * outlier_frac)
  outlier_idx <- if (n_out > 0) sample.int(N, n_out) else integer(0)
  is_outlier <- logical(N)
  is_outlier[outlier_idx] <- TRUE
  
  severity <- if (is.null(outlier_severity)) {
    switch(outlier_type, spike = 0.9, random_composition = 0.5, shock = 6)
  } else outlier_severity
  
  for (i in outlier_idx) {
    p <- switch(outlier_type,
                spike = {
                  p0 <- rep((1 - severity) / (Q - 1), Q)
                  p0[sample.int(Q, 1)] <- severity
                  p0
                },
                random_composition = {
                  g <- rgamma(Q, shape = severity, rate = 1)
                  g / sum(g)
                },
                shock = {
                  eta_i <- eta[i, ]
                  idx <- sample.int(Q, max(1, round(Q * 0.15)))
                  eta_i[idx] <- eta_i[idx] + rnorm(length(idx), 0, severity)
                  p_raw <- exp(eta_i - max(eta_i))
                  p_raw / sum(p_raw)
                })
    Y[i, ] <- as.vector(rmultinom(1, M[i], p))
  }
  
  list(Y = Y, X = X, group = group, M = M, N = N, Q = Q, P = P, J = J, K = K,
       n_per_group = n_j, is_outlier = is_outlier,
       true = list(mu = mu, phi = phi, B = B, F = Fm,
                   lambda = lambda, sigma2 = sigma2))
}