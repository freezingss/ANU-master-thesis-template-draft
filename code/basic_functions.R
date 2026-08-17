# Loading Matrix B
subspace_dist <- function(B1, B2) {
  q1 <- qr.Q(qr(B1)); q2 <- qr.Q(qr(B2))
  sv <- svd(crossprod(q1, q2))$d
  sv <- pmin(pmax(sv, 0), 1)
  sqrt(sum(1 - sv^2))
}

# Build Long Dataframe 
build_long <- function(Y, X, group) {
  N <- nrow(Y); Q <- ncol(Y); P <- ncol(X)
  long <- data.frame(
    count    = as.vector(Y),                                
    category = factor(rep(seq_len(Q), each = N)),
    obs      = factor(rep(seq_len(N), times = Q)),
    group    = factor(rep(group,      times = Q))
  )
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