source("basic_functions.R")

iso_newton_group <- function(Y_j, M_j, Eta0_j, tau2, lambda_init = NULL,
                             max_iter = 100, gtol = 1e-8, c1 = 1e-4) {
  Q <- ncol(Y_j)
  Nj <- nrow(Y_j)
  lambda <- if (is.null(lambda_init)) rep(0, Q) else lambda_init
  y_col <- colSums(Y_j)
  obj_fn <- function(lam) {
    Eta <- sweep(Eta0_j, 2, lam, "+")
    sum(Y_j * Eta) - sum(M_j * row_logsumexp(Eta)) - sum(lam^2) / (2 * tau2)
  }
  for (iter in 1:max_iter) {
    Eta <- sweep(Eta0_j, 2, lambda, "+")
    P <- t(row_softmax(Eta))
    Mp <- as.numeric(P %*% M_j)
    g <- y_col - Mp - lambda / tau2
    if (max(abs(g)) < gtol) break
    Dvec <- Mp + 1 / tau2
    Dinv_g <- g / Dvec
    Dinv_P <- P / Dvec
    Z <- diag(1 / M_j, Nj) - crossprod(P, Dinv_P)
    rhs <- crossprod(P, Dinv_g)
    corr <- tryCatch(solve(Z, rhs), error = function(e) rhs * 0)
    dir <- Dinv_g + as.numeric(Dinv_P %*% corr)
    gd <- sum(g * dir)
    obj0 <- obj_fn(lambda)
    step <- 1
    accepted <- FALSE
    repeat {
      if (obj_fn(lambda + step * dir) >= obj0 + c1 * step * gd) { accepted <- TRUE; break }
      step <- step / 2
      if (step < 1e-10) break
    }
    if (!accepted) {
      warning(sprintf("iso_newton_group: line search failed at inner iter %d (grad_norm=%.4g); returning current lambda", iter, max(abs(g))))
      break
    }
    lambda <- lambda + step * dir
  }
  Eta <- sweep(Eta0_j, 2, lambda, "+")
  P <- t(row_softmax(Eta))
  dF <- as.numeric(P %*% M_j)
  b <- (dF + 1 / tau2) * lambda - as.numeric(P %*% (M_j * as.numeric(crossprod(P, lambda))))
  list(lambda_hat = lambda, P = P, Mj = M_j, dF = dF, b = b, nj = Nj, tau2 = tau2)
}

estep_iso_wb <- function(Y, X, group, M, mu, phi, tau2, J = max(group),
                         max_iter = 100, gtol = 1e-8, lambda_warm = NULL) {
  Q <- ncol(Y)
  Xphi <- X %*% phi
  parts <- vector("list", J)
  lambda_hat <- matrix(0, Q, J)
  for (j in seq_len(J)) {
    idx <- which(group == j)
    Eta0 <- matrix(mu, length(idx), Q, byrow = TRUE) + Xphi[idx, , drop = FALSE]
    lin <- if (is.null(lambda_warm)) NULL else lambda_warm[, j]
    o <- iso_newton_group(Y[idx, , drop = FALSE], M[idx], Eta0, tau2,
                          lambda_init = lin, max_iter = max_iter, gtol = gtol)
    parts[[j]] <- o
    lambda_hat[, j] <- o$lambda_hat
  }
  list(lambda_hat = lambda_hat, parts = parts)
}

fhem_wb_parts <- function(Y, X, group, M, mu, phi, lambda_hat, tau2, J = max(group)) {
  Q <- ncol(Y)
  Xphi <- X %*% phi
  parts <- vector("list", J)
  for (j in seq_len(J)) {
    idx <- which(group == j)
    nj <- length(idx)
    lam <- lambda_hat[, j]
    Eta <- matrix(mu, nj, Q, byrow = TRUE) + Xphi[idx, , drop = FALSE] +
      matrix(lam, nj, Q, byrow = TRUE)
    Pj <- t(row_softmax(Eta))
    Mj <- M[idx]
    dF <- as.numeric(Pj %*% Mj)
    b <- (dF + 1 / tau2) * lam - as.numeric(Pj %*% (Mj * as.numeric(crossprod(Pj, lam))))
    parts[[j]] <- list(lambda_hat = lam, P = Pj, Mj = Mj, dF = dF, b = b,
                       nj = nj, tau2 = tau2)
  }
  parts
}

safe_chol <- function(A, jitter0 = 1e-10) {
  R <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(R)) return(R)
  s <- mean(abs(diag(A))) + 1e-300
  jit <- jitter0 * s
  for (k in 1:8) {
    R <- tryCatch(chol(A + diag(jit, nrow(A))), error = function(e) NULL)
    if (!is.null(R)) {
      warning(sprintf("safe_chol: matrix not numerically PD; added jitter %.3g", jit))
      return(R)
    }
    jit <- jit * 10
  }
  stop("safe_chol: Cholesky failed after jittering")
}

fhem_wb_const <- function(parts, ridge) {
  if (ridge <= 0)
    return(list(lt = vector("list", length(parts)), const = 0, exact = FALSE))
  lt <- vector("list", length(parts))
  tot <- 0
  for (j in seq_along(parts)) {
    pt <- parts[[j]]
    Dr <- pt$dF + ridge
    A <- pt$P / Dr
    Zf <- crossprod(pt$P, A)
    Zf <- diag(1 / pt$Mj, pt$nj) - (Zf + t(Zf)) / 2
    Rf <- safe_chol(Zf)
    logdetF <- sum(log(Dr)) + sum(log(pt$Mj)) + 2 * sum(log(diag(Rf)))
    v <- pt$lambda_hat / Dr +
      as.numeric(A %*% (chol2inv(Rf) %*% crossprod(A, pt$lambda_hat)))
    lt[[j]] <- pt$lambda_hat + (1 / pt$tau2 - ridge) * v
    tot <- tot + 0.5 * logdetF - 0.5 * sum(lt[[j]] * pt$b)
  }
  list(lt = lt, const = tot, exact = TRUE)
}

fhem_wb_G2 <- function(B, sigma2) sigma2^2 * diag(ncol(B)) + sigma2 * crossprod(B)

fhem_wb_estep <- function(parts, B, sigma2, ridge) {
  J <- length(parts)
  Q <- length(parts[[1]]$dF)
  K <- ncol(B)
  nj <- vapply(parts, function(p) p$nj, integer(1))
  wid <- nj + K
  off <- c(0, cumsum(wid))
  Tot <- off[J + 1]
  G2 <- fhem_wb_G2(B, sigma2)
  logdet_G2 <- 2 * sum(log(diag(safe_chol((G2 + t(G2)) / 2))))
  m <- matrix(0, Q, J)
  Astack <- matrix(0, Q, Tot)
  Cstack <- matrix(0, Q, Tot)
  Dsum <- numeric(Q)
  sum_log_Dt <- 0; sum_log_Mj <- 0; sum_logdetZ <- 0; quad <- 0
  for (j in seq_len(J)) {
    pt <- parts[[j]]
    n <- pt$nj
    Dt <- pt$dF + ridge + 1 / sigma2
    U <- cbind(pt$P, B)
    A <- U / Dt
    Z <- crossprod(U, A)
    Z <- -(Z + t(Z)) / 2
    Z[1:n, 1:n] <- Z[1:n, 1:n] + diag(1 / pt$Mj, n)
    Z[(n + 1):(n + K), (n + 1):(n + K)] <- Z[(n + 1):(n + K), (n + 1):(n + K)] + G2
    R <- safe_chol(Z)
    Zinv <- chol2inv(R)
    cols <- (off[j] + 1):off[j + 1]
    Astack[, cols] <- A
    Cstack[, cols] <- A %*% Zinv
    Dsum <- Dsum + 1 / Dt
    mj <- pt$b / Dt + as.numeric(Cstack[, cols, drop = FALSE] %*% crossprod(A, pt$b))
    m[, j] <- mj
    sum_log_Dt <- sum_log_Dt + sum(log(Dt))
    sum_log_Mj <- sum_log_Mj + sum(log(pt$Mj))
    sum_logdetZ <- sum_logdetZ + 2 * sum(log(diag(R)))
    quad <- quad + sum(pt$b * mj)
  }
  list(m = m, Astack = Astack, Cstack = Cstack, Dsum = Dsum, J = J, Q = Q, K = K,
       logdet_G2 = logdet_G2, sum_log_Dt = sum_log_Dt, sum_log_Mj = sum_log_Mj,
       sum_logdetZ = sum_logdetZ, quad = quad, off = off)
}

fhem_wb_logdet_Sigma <- function(sigma2, K, Q, logdet_G2) {
  Q * log(sigma2) + logdet_G2 - 2 * K * log(sigma2)
}

fhem_wb_loglik <- function(es, sigma2, const_tot) {
  ldS <- fhem_wb_logdet_Sigma(sigma2, es$K, es$Q, es$logdet_G2)
  ldV <- -es$sum_log_Dt - es$sum_log_Mj + es$J * es$logdet_G2 - es$sum_logdetZ
  -0.5 * es$J * ldS + 0.5 * ldV + 0.5 * es$quad + const_tot
}

fhem_wb_trS <- function(es, deflate = FALSE) {
  tot <- (sum(es$m^2) + sum(es$Dsum) + sum(es$Astack * es$Cstack)) / es$J
  if (!deflate) return(tot)
  cA <- colSums(es$Astack); cC <- colSums(es$Cstack)
  uSu <- (sum(colSums(es$m)^2) + sum(es$Dsum) + sum(cA * cC)) / (es$J * es$Q)
  tot - uSu
}

fhem_wb_Sop <- function(es, deflate = FALSE) {
  m <- es$m; A <- es$Astack; C <- es$Cstack; D <- es$Dsum; J <- es$J
  function(x, args = NULL) {
    if (deflate) x <- x - mean(x)
    v <- (as.numeric(m %*% crossprod(m, x)) + D * x +
            as.numeric(C %*% crossprod(A, x))) / J
    if (deflate) v <- v - mean(v)
    v
  }
}

fhem_wb_Sdense <- function(es, deflate = FALSE) {
  S <- (tcrossprod(es$m) + diag(es$Dsum) + tcrossprod(es$Cstack, es$Astack)) / es$J
  S <- (S + t(S)) / 2
  if (deflate) {
    r <- rowMeans(S)
    S <- S - outer(r, rep(1, es$Q)) - outer(rep(1, es$Q), r) + mean(S)
    S <- (S + t(S)) / 2
  }
  S
}

fhem_lanczos_topk <- function(op, n, k, m = NULL, reltol = 1e-11) {
  if (is.null(m)) m <- min(n, max(6L * k + 40L, 60L))
  v0 <- sin(seq_len(n) * 0.7853981633974483) + cos(seq_len(n) * 0.3183098861837907)
  run <- function(m) {
    V <- matrix(0, n, m)
    alpha <- numeric(m); beta <- numeric(m)
    V[, 1] <- v0 / sqrt(sum(v0^2))
    mm <- m
    for (i in seq_len(m)) {
      w <- op(V[, i])
      alpha[i] <- sum(w * V[, i])
      Vi <- V[, 1:i, drop = FALSE]
      w <- w - Vi %*% crossprod(Vi, w)
      w <- w - Vi %*% crossprod(Vi, w)
      bn <- sqrt(sum(w^2))
      if (i == m || bn < 1e-13) { mm <- i; break }
      beta[i] <- bn
      V[, i + 1] <- w / bn
    }
    Tm <- diag(alpha[1:mm], mm)
    if (mm > 1) {
      id <- cbind(2:mm, 1:(mm - 1))
      Tm[id] <- beta[1:(mm - 1)]
      Tm[id[, c(2, 1)]] <- beta[1:(mm - 1)]
    }
    eg <- eigen(Tm, symmetric = TRUE)
    ord <- order(eg$values, decreasing = TRUE)[1:k]
    list(values = eg$values[ord],
         vectors = V[, 1:mm, drop = FALSE] %*% eg$vectors[, ord, drop = FALSE],
         resid = if (mm < m) rep(0, k) else abs(beta[mm - 1] * eg$vectors[mm, ord]))
  }
  out <- run(m)
  tries <- 0
  while (max(out$resid) > reltol * max(abs(out$values)) && m < n && tries < 3) {
    m <- min(n, 2L * m)
    out <- run(m)
    tries <- tries + 1
  }
  out
}

fhem_wb_topk <- function(es, K, method, deflate = FALSE) {
  Q <- es$Q
  if (method == "auto") method <- if (Q <= 150) "dense" else "lanczos"
  if (method == "dense") {
    S <- fhem_wb_Sdense(es, deflate)
    eg <- eigen(S, symmetric = TRUE)
    return(list(values = eg$values[1:K], vectors = eg$vectors[, 1:K, drop = FALSE],
                trS = sum(diag(S))))
  }
  op <- fhem_wb_Sop(es, deflate)
  trS <- fhem_wb_trS(es, deflate)
  if (requireNamespace("RSpectra", quietly = TRUE) && K + 2 < Q) {
    eig <- tryCatch(RSpectra::eigs_sym(op, k = K, n = Q, which = "LA",
                                       opts = list(tol = 1e-14, maxitr = 5000)),
                    error = function(e) NULL)
    if (!is.null(eig) && eig$nconv >= K)
      return(list(values = eig$values, vectors = eig$vectors, trS = trS))
  }
  lz <- fhem_lanczos_topk(op, Q, K)
  list(values = lz$values, vectors = lz$vectors, trS = trS)
}

fhem_wb_mstep <- function(es, K, method, sigma2_floor = 1e-8, deflate = FALSE) {
  tk <- fhem_wb_topk(es, K, method, deflate)
  lam_K <- tk$values
  df <- if (deflate) es$Q - 1 - K else es$Q - K
  sigma2 <- if (df > 0) max((tk$trS - sum(lam_K)) / df, sigma2_floor) else sigma2_floor
  B <- tk$vectors %*% diag(sqrt(pmax(lam_K - sigma2, 0)), K)
  list(B = B, sigma2 = sigma2, trS = tk$trS, lam_K = lam_K)
}

fhem_gaussian_em_wb <- function(parts, const_tot, tau2, K, ridge,
                                max_iter = 300, tol = 1e-6,
                                method = "auto", sigma2_floor = 1e-8,
                                deflate = FALSE) {
  Q <- length(parts[[1]]$dF)
  B <- matrix(0, Q, K)
  sigma2 <- tau2
  es <- fhem_wb_estep(parts, B, sigma2, ridge)
  ll <- fhem_wb_loglik(es, sigma2, const_tot)
  trace <- ll
  it <- 0
  for (it in 1:max_iter) {
    ms <- fhem_wb_mstep(es, K, method, sigma2_floor, deflate)
    B <- ms$B
    sigma2 <- ms$sigma2
    es <- fhem_wb_estep(parts, B, sigma2, ridge)
    ll_new <- fhem_wb_loglik(es, sigma2, const_tot)
    trace <- c(trace, ll_new)
    if (abs(ll_new - ll) < tol) break
    ll <- ll_new
  }
  list(B = B, sigma2 = sigma2, es = es, trace = trace, iters = it,
       monotone = all(diff(trace) > -1e-6))
}

fhem_from_parts_wb <- function(parts, tau2, K, ridge = 1e-3,
                               max_iter = 300, tol = 1e-6,
                               method = "auto", sigma2_floor = 1e-8,
                               deflate = FALSE) {
  cn <- fhem_wb_const(parts, ridge)
  em <- fhem_gaussian_em_wb(parts, cn$const, tau2, K, ridge,
                            max_iter = max_iter, tol = tol, method = method,
                            sigma2_floor = sigma2_floor, deflate = deflate)
  list(B = em$B, sigma2 = em$sigma2, lt = cn$lt, m = em$es$m,
       ll_exact = cn$exact, ridge = ridge,
       trace = em$trace, iters = em$iters, monotone = em$monotone)
}

fhem_from_iso_estep_wb <- function(lambda_hat, Y, X, group, M, mu, phi, tau2, K,
                                   ridge = 1e-3, max_iter = 300, tol = 1e-6,
                                   method = "auto", sigma2_floor = 1e-8,
                                   deflate = FALSE) {
  parts <- fhem_wb_parts(Y, X, group, M, mu, phi, lambda_hat, tau2)
  fhem_from_parts_wb(parts, tau2, K, ridge = ridge, max_iter = max_iter, tol = tol,
                     method = method, sigma2_floor = sigma2_floor, deflate = deflate)
}

fit_pfa_fhem_wb <- function(Y, X, group, M, K,
                        mu0, phi0, B0, sigma2_0,
                        estep_iso_fn = NULL, mstep_phi_fn,
                        apply_PLT_fn = apply_PLT,
                        ridge = 1e-3, max_iter = 30, tol_outer = 1e-6,
                        fhem_max_iter = 300, fhem_tol = 1e-6,
                        estep_max_iter = 100, estep_gtol = 1e-8,
                        method = "auto", sigma2_floor = 1e-8, deflate = FALSE,
                        verbose = TRUE, trace = TRUE, B_true = NULL) {
  Q <- length(mu0)
  J <- max(group)
  mu <- mu0
  phi <- phi0
  B <- B0
  sigma2 <- sigma2_0
  tau2 <- sum(B^2) / Q + sigma2
  hist_df <- data.frame(iter = integer(0), tau2 = numeric(0), sigma2 = numeric(0),
                        fh_ll = numeric(0), inner_iters = integer(0))
  B_dist_trace <- numeric(0)
  u_frac_trace <- numeric(0)
  ll_prev <- -Inf
  converged <- FALSE
  best_iter <- NA_integer_
  best_ll <- -Inf
  best_state <- NULL

  get_parts <- function(mu, phi, tau2) {
    if (is.null(estep_iso_fn)) {
      e <- estep_iso_wb(Y, X, group, M, mu, phi, tau2, J,
                        max_iter = estep_max_iter, gtol = estep_gtol)
      list(lambda_hat = e$lambda_hat, parts = e$parts)
    } else {
      e <- estep_iso_fn(mu, phi, tau2)
      list(lambda_hat = e$lambda_hat,
           parts = fhem_wb_parts(Y, X, group, M, mu, phi, e$lambda_hat, tau2, J))
    }
  }

  for (outer in 1:max_iter) {
    es <- get_parts(mu, phi, tau2)
    mp <- mstep_phi_fn(es$lambda_hat, mu, phi)
    mu <- mp$mu; phi <- mp$phi
    fh <- fhem_from_parts_wb(es$parts, tau2, K, ridge = ridge,
                             max_iter = fhem_max_iter, tol = fhem_tol,
                             method = method, sigma2_floor = sigma2_floor,
                             deflate = deflate)
    B <- apply_PLT_fn(fh$B)
    sigma2 <- fh$sigma2
    tau2_new <- sum(B^2) / Q + sigma2
    ll <- tail(fh$trace, 1)
    hist_df <- rbind(hist_df, data.frame(iter = outer, tau2 = tau2, sigma2 = sigma2,
                                         fh_ll = ll, inner_iters = fh$iters))
    if (trace && !is.null(B_true)) B_dist_trace <- c(B_dist_trace, subspace_distance(B, B_true))
    if (trace) u_frac_trace <- c(u_frac_trace, u_frac(B))
    if (trace && ll > best_ll) {
      best_ll <- ll; best_iter <- outer
      best_state <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2, tau2 = tau2_new,
                         lambda_hat = es$lambda_hat, lt = fh$lt)
    }
    if (verbose) cat(sprintf(
      "outer %2d  tau2=%.4f  sigma2=%.4f  fh_ll=%.2f  inner_iters=%d  inner_monotone=%s\n",
      outer, tau2, sigma2, ll, fh$iters, fh$monotone))
    converged <- outer > 1 && abs(ll - ll_prev) < tol_outer
    ll_prev <- ll
    tau2 <- tau2_new
    if (converged) break
  }

  es_final <- get_parts(mu, phi, tau2)
  fh_final <- fhem_from_parts_wb(es_final$parts, tau2, K, ridge = ridge,
                                 max_iter = fhem_max_iter, tol = fhem_tol,
                                 method = method, sigma2_floor = sigma2_floor,
                                 deflate = deflate)
  B <- apply_PLT_fn(fh_final$B)
  sigma2 <- fh_final$sigma2
  tau2 <- sum(B^2) / Q + sigma2

  out <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2, tau2 = tau2,
              lt = fh_final$lt, m_hat = fh_final$m,
              lambda_hat = es_final$lambda_hat,
              trace = hist_df, fh_ll_final = tail(fh_final$trace, 1),
              fh_trace_final = fh_final$trace, fh_ll_exact = fh_final$ll_exact,
              converged = converged, iterations = nrow(hist_df),
              n_iter = nrow(hist_df), sigma2_trace = hist_df$sigma2,
              log_evidence = hist_df$fh_ll)
  if (trace && !is.null(B_true)) out$B_dist_trace <- B_dist_trace
  if (trace) out$u_frac_trace <- u_frac_trace
  if (trace) out <- c(out, list(best_iter = best_iter, best_ll = best_ll, best = best_state))
  out
}

fit_pfa_fhem <- function(Y, X, group, K,
                          M = rowSums(Y),
                          max_iter = 30, tol = 1e-6,
                          lambda_phi = 0, sigma2_init = 0.3,
                          verbose = FALSE,
                          estep_max_iter = 100, estep_gtol = 1e-8,
                          fhem_max_iter = 300, fhem_tol = 1e-6,
                          ridge = 1e-3,
                          method = c("auto", "lanczos", "dense"),
                          sigma2_floor = 1e-8, deflate = FALSE,
                          external_estep = FALSE,
                          B_true = NULL, trace = TRUE) {
  method <- match.arg(method)
  Q <- ncol(Y)
  J <- max(group)
  th <- init_theta(Y, X, group, K, sigma2_init)
  estep_iso_fn <- if (external_estep) {
    function(mu, phi, tau2)
      estep_wbonly(J, group, Y, X, M, mu, phi, matrix(0, Q, K), tau2,
                   estep_max_iter, estep_gtol)
  } else NULL
  mstep_phi_fn <- function(lambda_hat, mu, phi)
    mstep_phi_wbonly(Y, X, group, lambda_hat, mu, phi, lambda_phi = lambda_phi)

  fit_pfa_fhem_wb(Y, X, group, M, K,
                mu0 = th$mu, phi0 = th$phi, B0 = th$B, sigma2_0 = th$sigma2,
                estep_iso_fn = estep_iso_fn, mstep_phi_fn = mstep_phi_fn,
                apply_PLT_fn = apply_PLT,
                ridge = ridge, max_iter = max_iter, tol_outer = tol,
                fhem_max_iter = fhem_max_iter, fhem_tol = fhem_tol,
                estep_max_iter = estep_max_iter, estep_gtol = estep_gtol,
                method = method, sigma2_floor = sigma2_floor, deflate = deflate,
                verbose = verbose, trace = trace, B_true = B_true)
}
