source("setup.R")

# Long Table

build_long <- function(Y, X, group) {
  N <- nrow(Y)
  Q <- ncol(Y)
  P <- ncol(X)
  long <- data.frame(
    count     = as.vector(Y),
    category  = factor(rep(seq_len(Q), each = N)),
    group     = factor(rep(group, times = Q)),
    log_total = rep(log(rowSums(Y)), times = Q)
  )
  long$obs <- interaction(long$group, long$category, drop = TRUE)
  if (P > 1) {
    for (p in 2:P) long[[paste0("x", p)]] <- rep(X[, p], times = Q)
  }
  long
}

apply_PLT <- function(B) {
  # Add: remove the average
  B <- sweep(B, 2, colMeans(B))
  K <- ncol(B)
  qr_obj <- qr(t(B[1:K, , drop = FALSE]))
  B_rot <- B %*% qr.Q(qr_obj)
  for (k in 1:K) if (B_rot[k, k] < 0) B_rot[, k] <- -B_rot[, k]
  for (k in 1:K) if (k > 1) B_rot[1:(k - 1), k] <- 0
  B_rot
}

init_theta <- function(Y, X, group, K, sigma2_init) {
  Q <- ncol(Y); P <- ncol(X); J <- max(group)
  ap <- colMeans(Y / pmax(rowSums(Y), 1))
  mu <- log(ap + 1e-8); mu <- mu - mean(mu)
  gm <- matrix(0, J, Q)
  for (j in 1:J) {
    idx <- which(group == j)
    if (length(idx)) {
      gs <- colSums(Y[idx, , drop = FALSE])
      gm[j, ] <- log((gs + 1e-5) / (sum(gs) + Q * 1e-5))
    }
  }
  sv0 <- svd(t(sweep(gm, 2, colMeans(gm))), nu = K, nv = K)
  list(mu = mu, phi = matrix(0, P, Q),
       B = apply_PLT(sv0$u %*% diag(pmax(sv0$d[1:K] * 0.5, 0.1), K)),
       sigma2 = sigma2_init)
}

lambda_mean_cor <- function(lambda_hat, lambda_true) {
  J <- ncol(lambda_true)
  mean(sapply(1:J, 
              function(j) cor(lambda_hat[, j], 
                              lambda_true[, j])))
}

all_perms <- function(K) {
  if (K == 1) return(matrix(1L, 1, 1))
  prev <- all_perms(K - 1)
  out <- matrix(NA_integer_, nrow(prev) * K, K)
  row <- 1
  for (i in 1:nrow(prev)) {
    for (pos in 1:K) {
      out[row, ] <- append(prev[i, ], K, after = pos - 1)
      row <- row + 1
    }
  }
  out
}

# Metrices: Determine the Accuracy
subspace_distance <- function(B_hat, B_true, normalize = FALSE) {
  K <- ncol(B_true)
  Q1 <- qr.Q(qr(B_hat))
  Q2 <- qr.Q(qr(B_true))
  sv <- svd(crossprod(Q1, Q2))$d
  d <- sqrt(2 * sum(1 - pmin(sv, 1)^2))
  if (normalize) d / sqrt(2 * K) else d
}

principal_angles <- function(B1, B2, degrees = TRUE) {
  K <- ncol(B1)
  U1 <- qr.Q(qr(B1))[, 1:K, drop = FALSE]
  U2 <- qr.Q(qr(B2))[, 1:K, drop = FALSE]
  s <- pmin(pmax(svd(crossprod(U1, U2))$d, -1), 1)
  th <- acos(s)
  if (degrees) th * 180 / pi else th
}

procrustes_align <- function(B_hat, B_true, allow_scale = FALSE) {
  sv <- svd(crossprod(B_hat, B_true))
  Bal <- B_hat %*% (sv$u %*% t(sv$v))
  if (allow_scale)
    Bal <- Bal * (sum(Bal * B_true) / sum(Bal * Bal))
  Bal
}

tucker_congruence <- function(B_hat, B_true, align = TRUE) {
  Bal <- if (align) procrustes_align(B_hat, B_true) else B_hat
  K <- ncol(B_true)
  phi <- sapply(1:K, function(k) {
    a <- Bal[, k]
    b <- B_true[, k]
    d <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (d < 1e-12) NA_real_ else sum(a * b) / d
  })
  list(per_factor = phi, mean = mean(abs(phi), na.rm = TRUE))
}

rv_coefficient <- function(B_hat, B_true) {
  A <- tcrossprod(B_hat)
  Bm <- tcrossprod(B_true)
  sum(A * Bm) / sqrt(sum(A * A) * sum(Bm * Bm))
}

bbp_floor <- function(B_true, sigma2_true, J, normalize = FALSE) {
  Q <- nrow(B_true)
  K <- ncol(B_true)
  gamma <- Q / J
  theta <- eigen(tcrossprod(B_true), symmetric = TRUE)$values[1:K] / sigma2_true
  d_inf <- sqrt(2 * sum(sapply(theta, function(th) {
    if (th <= sqrt(gamma)) return(1)
    1 - (1 - gamma / th^2) / (1 + gamma / th)
  })))
  if (normalize) d_inf / sqrt(2 * K) else d_inf
}

loading_rmse <- function(B_hat, B_true, align = TRUE, allow_scale = FALSE) {
  D <- (if (align) procrustes_align(B_hat, B_true, allow_scale) else B_hat) - B_true
  sqrt(mean(D^2))
}

sigma_error <- function(B_hat, sigma2_hat, B_true, sigma2_true, relative = TRUE) {
  Sh <- tcrossprod(B_hat) + sigma2_hat * diag(nrow(B_hat))
  St <- tcrossprod(B_true) + sigma2_true * diag(nrow(B_true))
  num <- sqrt(sum((Sh - St)^2))
  if (relative) num / sqrt(sum(St^2)) else num
}

B_metrics <- function(B, sigma2, B_true, sigma2_true) {
  if (is.null(B))
    return(c(d_true = NA, ang_max = NA, tucker = NA, rv = NA, rmse = NA, sig_err = NA))
  c(d_true  = subspace_distance(B, B_true),
    ang_max = max(principal_angles(B, B_true)),
    tucker = tucker_congruence(B, B_true)$mean,
    rv = rv_coefficient(B, B_true),
    rmse = loading_rmse(B, B_true),
    sig_err = sigma_error(B, sigma2, B_true, sigma2_true))
}

row_softmax <- function(eta) {
  mx <- apply(eta, 1, max)
  ee <- exp(eta - mx)
  ee / rowSums(ee)
}
row_logsumexp <- function(eta) {
  mx <- apply(eta, 1, max)
  mx + log(rowSums(exp(eta - mx)))
}

# # Classic PPCA, no deflate
# ppca_closed <- function(S, K) {
#   Q <- nrow(S)
#   e <- eigen((S + t(S)) / 2, symmetric = TRUE)
#   ev <- e$values
#   q <- length(ev)
#   sigma2 <- if (K < q) max(mean(ev[(K + 1):q]), 1e-8) else 1e-8
#   B <- e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(ev[1:K] - sigma2, 0)), K)
#   list(B = B, sigma2 = sigma2, ev = ev, Sigma = tcrossprod(B) + sigma2 * diag(q))
# }

# Adjusted PPCA on Q-1 dimensions
ppca_closed <- function(S, K) {
  Q <- nrow(S)
  V <- qr.Q(qr(matrix(1, Q, 1)), complete = TRUE)[, 2:Q, drop = FALSE]

  S_perp <- t(V) %*% S %*% V
  S_perp <- (S_perp + t(S_perp)) / 2

  e <- eigen(S_perp, symmetric = TRUE)
  ev <- e$values
  q <- length(ev)
  sigma2 <- if (K < q) max(mean(ev[(K + 1):q]), 1e-8) else 1e-8
  B_perp <- e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(ev[1:K] - sigma2, 0)), K)

  B <- V %*% B_perp

  list(B = B, sigma2 = sigma2, ev = ev, Sigma = tcrossprod(B) + sigma2 * diag(Q))
}

correct_lambda_edgeworth <- function(Q, J, lambda_hat, S_hat, Y, X, group, M, mu, phi) {
  lam_corr <- lambda_hat
  rel_size <- numeric(J)
  Xphi <- X %*% phi
  for (j in seq_len(J)) {
    idx <- which(group == j)
    n_j <- length(idx)
    Sj <- S_hat[[j]]
    diagS <- diag(Sj)
    
    eta <- matrix(mu, n_j, Q, byrow = TRUE) +
      Xphi[idx, , drop = FALSE] +
      matrix(lambda_hat[, j], n_j, Q, byrow = TRUE)
    eta <- eta - apply(eta, 1, max)
    pi_mat <- exp(eta)
    pi_mat <- pi_mat / rowSums(pi_mat)
    
    TS <- numeric(Q)
    for (i in seq_len(n_j)) {
      pi_i <- pi_mat[i, ]
      Spi_i <- as.vector(Sj %*% pi_i)
      piDiagS <- sum(pi_i * diagS)
      piSpi <- sum(pi_i * Spi_i)
      TS <- TS + M[idx[i]] * pi_i * (diagS - 2 * Spi_i - piDiagS + 2 * piSpi)
    }
    
    mu1 <- -0.5 * as.vector(Sj %*% TS)
    
    lam_corr[, j] <- lambda_hat[, j] + mu1
    rel_size[j] <- sqrt(sum(mu1^2)) / max(sqrt(sum(lambda_hat[, j]^2)), 1e-8)
  }
  
  list(lambda_corrected = lam_corr, rel_size = rel_size)
}

apply_lambda_correction <- function(Q, J, lambda_hat, S_hat, Y, X, group, M, mu, phi, corr_max_rel) {
  cc <- correct_lambda_edgeworth(Q, J, lambda_hat, S_hat, Y, X, group, M, mu, phi)
  lam_c <- cc$lambda_corrected
  n_capped <- 0L
  for (j in seq_len(J)) {
    if (is.finite(cc$rel_size[j]) && cc$rel_size[j] > corr_max_rel) {
      mu1 <- lam_c[, j] - lambda_hat[, j]
      lam_c[, j] <- lambda_hat[, j] + mu1 * (corr_max_rel / cc$rel_size[j])
      n_capped <- n_capped + 1L
    }
  }
  list(lambda_corrected = lam_c, rel_size = cc$rel_size, n_capped = n_capped)
}