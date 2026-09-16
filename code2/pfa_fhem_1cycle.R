source("basic_functions.R")

iso_newton_group <- function(Y_j, M_j, Eta0_j, tau2, lambda_init = NULL, max_iter = 50, gtol = 1e-8, c1 = 1e-4) {
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
      if (step < 1e-6) break
    }
    if (!accepted) {
      warning(sprintf("iso_newton_group: line search failed to find an accepted step at inner iter %d (grad_norm=%.4g); returning current unconverged lambda", iter, max(abs(g))))
      break
    }
    lambda <- lambda + step * dir
  }
  Eta <- sweep(Eta0_j, 2, lambda, "+")
  P <- t(row_softmax(Eta))
  Mp <- as.numeric(P %*% M_j)
  D_F <- Mp
  D_iso <- D_F + 1 / tau2
  b <- D_iso * lambda - as.numeric(P %*% (M_j * as.numeric(crossprod(P, lambda))))
  list(lambda_hat = lambda, P = P, M_j = M_j, D_F = D_F, b = b)
}

fhem_estep_group <- function(iso_out, B, sigma2) {
  Q <- length(iso_out$D_F)
  K <- ncol(B)
  Nj <- length(iso_out$M_j)
  stopifnot(nrow(B) == Q)
  Dtilde <- iso_out$D_F + 1 / sigma2
  Utilde <- cbind(iso_out$P, B)
  BtB <- crossprod(B)
  Wtilde_inv <- matrix(0, Nj + K, Nj + K)
  Wtilde_inv[1:Nj, 1:Nj] <- diag(1 / iso_out$M_j, Nj)
  Wtilde_inv[(Nj + 1):(Nj + K), (Nj + 1):(Nj + K)] <- sigma2^2 * diag(K) + sigma2 * BtB
  Z <- Wtilde_inv - crossprod(Utilde, Utilde / Dtilde)
  Zinv <- tryCatch(solve(Z), error = function(e) {
    warning(sprintf("fhem_estep_group: solve(Z) failed (%s); falling back to a heavily damped approximation", conditionMessage(e)))
    solve(Z + 1e-3 * diag(nrow(Z)))
  })
  Dinv_b <- iso_out$b / Dtilde
  m <- Dinv_b + (Utilde %*% (Zinv %*% crossprod(Utilde, Dinv_b))) / Dtilde
  list(Dtilde = Dtilde, Utilde = Utilde, Zinv = Zinv, m = as.vector(m))
}

apply_Vj_corr <- function(Dtilde, Utilde, Zinv, x) {
  Dinv_x <- x / Dtilde
  corr <- Zinv %*% crossprod(Utilde, Dinv_x)
  (Utilde %*% corr) / Dtilde
}

build_S_operator <- function(estep_list, J) {
  Mmat <- sapply(estep_list, function(e) e$m)
  Dbar <- Reduce(`+`, lapply(estep_list, function(e) 1 / e$Dtilde)) / J
  function(x, args) {
    term1 <- Mmat %*% crossprod(Mmat, x) / J
    term2 <- Dbar * x
    term3 <- Reduce(`+`, lapply(estep_list, function(e) apply_Vj_corr(e$Dtilde, e$Utilde, e$Zinv, x))) / J
    as.numeric(term1 + term2 + term3)
  }
}

fhem_trS <- function(estep_list, J) {
  mean(sapply(estep_list, function(e) {
    sum(e$m^2) + sum(1 / e$Dtilde) + sum(diag(e$Zinv %*% crossprod(e$Utilde / e$Dtilde)))
  }))
}

Vj_dense <- function(Dtilde, Utilde, Zinv) {
  Dinv <- 1 / Dtilde
  diag(Dinv) + (Dinv * Utilde) %*% Zinv %*% t(Dinv * Utilde)
}

fhem_mstep <- function(estep_list, J, Q, K, sigma2_floor = 1e-8, method = c("lanczos", "dense")) {
  method <- match.arg(method)
  trS <- fhem_trS(estep_list, J)
  if (method == "lanczos") {
    Sv_op <- build_S_operator(estep_list, J)
    eig <- RSpectra::eigs_sym(Sv_op, k = K, n = Q, which = "LA")
    if (eig$nconv < K)
      warning(sprintf("fhem_mstep: RSpectra::eigs_sym only converged %d of %d requested eigenvalues", eig$nconv, K))
    lam_K <- eig$values
    U_K <- eig$vectors
  } else {
    Mmat <- sapply(estep_list, function(e) e$m)
    S_dense <- Mmat %*% t(Mmat) / J +
      Reduce(`+`, lapply(estep_list, function(e) Vj_dense(e$Dtilde, e$Utilde, e$Zinv))) / J
    eg <- eigen(S_dense, symmetric = TRUE)
    lam_K <- eg$values[1:K]
    U_K <- eg$vectors[, 1:K, drop = FALSE]
  }
  sigma2_new <- if (K < Q) max((trS - sum(lam_K)) / (Q - K), sigma2_floor) else sigma2_floor
  B_new <- U_K %*% diag(sqrt(pmax(lam_K - sigma2_new, 0)), K)
  B_new <- apply_PLT(B_new)
  list(B = B_new, sigma2 = sigma2_new, trS = trS, lam_K = lam_K)
}

fit_fhem_woodbury <- function(Y, X, group, J, mu, phi, K, B_init, sigma2_init,
                              max_iter = 100, iso_gtol = 1e-8, estep_max_iter = 50,
                              sigma2_floor = 1e-8,
                              tau2_min = 1e-6, tol = 1e-8, method = c("lanczos", "dense"),
                              lambda_phi = 0, B_true = NULL, trace = TRUE) {
  method <- match.arg(method)
  Q <- ncol(Y)
  B <- B_init
  sigma2 <- sigma2_init
  sigma2_trace <- numeric(max_iter)
  B_dist_trace <- if (trace && !is.null(B_true)) numeric(max_iter) else NULL
  converged <- FALSE
  it <- 0

  for (it in 1:max_iter) {
    tau2 <- max(sum(B^2) / Q + sigma2, tau2_min)
    iso_list <- vector("list", J)
    for (j in 1:J) {
      idx <- which(group == j)
      Y_j <- Y[idx, , drop = FALSE]
      M_j <- rowSums(Y_j)
      Eta0_j <- matrix(mu, length(idx), Q, byrow = TRUE) + X[idx, , drop = FALSE] %*% phi
      iso_list[[j]] <- iso_newton_group(Y_j, M_j, Eta0_j, tau2, max_iter = estep_max_iter, gtol = iso_gtol)
    }
    estep_list <- vector("list", J)
    for (j in 1:J) {
      estep_list[[j]] <- fhem_estep_group(iso_list[[j]], B, sigma2)
    }

    lambda_hat_m <- sapply(estep_list, function(e) e$m)
    mp <- mstep_phi_wbonly(Y, X, group, lambda_hat_m, mu, phi, lambda_phi = lambda_phi)
    mu <- mp$mu
    phi <- mp$phi

    m_out <- fhem_mstep(estep_list, J, Q, K, sigma2_floor, method)
    sigma2_trace[it] <- m_out$sigma2
    delta_sigma2 <- abs(m_out$sigma2 - sigma2)
    delta_B <- max(abs(m_out$B - B))
    B <- m_out$B
    sigma2 <- m_out$sigma2
    if (trace && !is.null(B_true)) B_dist_trace[it] <- subspace_distance(B, B_true)
    if (delta_sigma2 < tol && delta_B < 1e-6) { converged <- TRUE; break }
  }

  tau2_final <- max(sum(B^2) / Q + sigma2, tau2_min)
  iso_list_final <- vector("list", J)
  for (j in 1:J) {
    idx <- which(group == j)
    Y_j <- Y[idx, , drop = FALSE]
    M_j <- rowSums(Y_j)
    Eta0_j <- matrix(mu, length(idx), Q, byrow = TRUE) + X[idx, , drop = FALSE] %*% phi
    iso_list_final[[j]] <- iso_newton_group(Y_j, M_j, Eta0_j, tau2_final, max_iter = estep_max_iter, gtol = iso_gtol)
  }
  estep_list_final <- vector("list", J)
  for (j in 1:J) estep_list_final[[j]] <- fhem_estep_group(iso_list_final[[j]], B, sigma2)
  lambda_hat <- sapply(iso_list_final, function(x) x$lambda_hat)
  m_hat <- sapply(estep_list_final, function(x) x$m)

  out <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2, tau2 = tau2_final,
             converged = converged,
             sigma2_trace = sigma2_trace[1:it], n_iter = it,
             lambda_hat = lambda_hat, m_hat = m_hat)
  if (trace && !is.null(B_true)) out$B_dist_trace <- B_dist_trace[1:it]
  out
}

fit_pfa_fh_em <- function(Y, X, group, K,
                          M = rowSums(Y),
                          max_iter = 100, tol = 1e-4,
                          lambda_phi = 0, sigma2_init = 0.3,
                          verbose = FALSE,
                          estep_max_iter = 100, estep_gtol = 1e-3,
                          method = c("lanczos", "dense"),
                          sigma2_floor = 1e-8, tau2_min = 1e-6,
                          B_true = NULL, trace = TRUE) {
  method <- match.arg(method)
  J <- max(group)
  th <- init_theta(Y, X, group, K, sigma2_init)
  fit <- fit_fhem_woodbury(Y, X, group, J, th$mu, th$phi, K, th$B, th$sigma2,
                           max_iter = max_iter, iso_gtol = estep_gtol,
                           estep_max_iter = estep_max_iter,
                           sigma2_floor = sigma2_floor, tau2_min = tau2_min,
                           tol = tol, method = method, lambda_phi = lambda_phi,
                           B_true = B_true, trace = trace)
  if (verbose) {
    bd <- if (!is.null(fit$B_dist_trace)) sprintf("  final B_dist=%.4f", tail(fit$B_dist_trace, 1)) else ""
    cat(sprintf("fit_pfa_fh_em: n_iter=%d  converged=%s  sigma2=%.4f%s\n",
                fit$n_iter, fit$converged, fit$sigma2, bd))
  }
  fit$iterations <- fit$n_iter
  fit
}
