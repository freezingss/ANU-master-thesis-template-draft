# Loading Matrix B
subspace_dist <- function(B1, B2) {
  q1 <- qr.Q(qr(B1)); q2 <- qr.Q(qr(B2))
  sv <- svd(crossprod(q1, q2))$d
  sv <- pmin(pmax(sv, 0), 1)
  sqrt(sum(1 - sv^2))
}

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

# Lower Triangular Transformation
apply_PLT <- function(B) {
  K <- ncol(B)
  qr_obj <- qr(t(B[1:K, , drop = FALSE]))
  B_rot <- B %*% qr.Q(qr_obj)
  for (k in 1:K) if (B_rot[k, k] < 0) B_rot[, k] <- -B_rot[, k]
  for (k in 1:K) if (k > 1) B_rot[1:(k - 1), k] <- 0
  B_rot
}

# Vector: Determine the Accuracy
lambda_mean_cor <- function(lambda_hat, lambda_true) {
  J <- ncol(lambda_true)
  mean(sapply(1:J, function(j) cor(lambda_hat[, j], lambda_true[, j])))
}

# Matrices: Determine the Accuracy
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

# Softmax
row_softmax <- function(eta) { 
  mx <- apply(eta, 1, max)
  ee <- exp(eta - mx)
  ee / rowSums(ee) 
}
row_logsumexp <- function(eta) { 
  mx <- apply(eta, 1, max)
  mx + log(rowSums(exp(eta - mx))) 
}

# Table: Determine converged or not
summarize_fit_traced <- function(fit, Q, K, J, N_per_group, seed, max_iter, tol) {
  data.frame(
    Q = Q, K = K, J = J, N_per_group = N_per_group, seed = seed,
    max_iter = max_iter, tol = tol,
    converged = fit$converged,
    iterations = fit$iterations,
    log_ev_final = tail(fit$log_evidence, 1),
    sigma2_final = fit$sigma2,
    # B_dist_final = if (!is.null(fit$B_dist_trace)) tail(fit$B_dist_trace, 1) else NA
    # The correlation between estimated B and true B
    
  )
}

#
