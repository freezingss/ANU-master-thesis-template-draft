simulate_pfa_data <- function(Q = 50, K = 2, J = 30, N_per_group = 15,
                              P = 3, sigma2 = 0.3, M_rate = 150,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  mu <- rnorm(Q, 0, 0.4); mu <- mu - mean(mu)
  phi <- matrix(0, P, Q)
  if (P > 1) phi[2:P, ] <- rnorm((P - 1) * Q, 0, 0.3)

  B <- apply_PLT(matrix(rnorm(Q * K, 0, 0.6), Q, K))
  Fm <- matrix(rnorm(J * K), J, K)                    
  lambda <- B %*% t(Fm) + matrix(rnorm(Q * J, 0, sqrt(sigma2)), Q, J)   

  N <- J * N_per_group
  group <- rep(1:J, each = N_per_group)
  X <- cbind(1, matrix(rnorm(N * (P - 1)), N, P - 1))
  M <- 1 + rpois(N, M_rate)

  eta <- sweep(X %*% phi, 2, mu, "+") + t(lambda[, group, drop = FALSE])
  Y <- t(sapply(1:N, function(i) {
    p <- exp(eta[i, ] - max(eta[i, ]))
    as.vector(rmultinom(1, M[i], p / sum(p)))
  }))

  list(Y = Y, X = X, group = group, M = M, N = N, Q = Q, P = P, J = J, K = K,
       true = list(mu = mu, phi = phi, B = B, F = Fm,
                   lambda = lambda, sigma2 = sigma2))
}



lambda_mean_cor <- function(lambda_hat, lambda_true) {
  J <- ncol(lambda_true)
  mean(sapply(1:J, function(j) cor(lambda_hat[, j], lambda_true[, j])))
}
