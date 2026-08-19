fit_pfa_woodbury_traced <- function(Y, X, group, K,
                                    M = rowSums(Y),
                                    max_iter = 60, tol = 1e-4,
                                    lambda_phi = 0, sigma2_init = 0.3,
                                    verbose = FALSE,
                                    estep_max_iter = 100, estep_gtol = 1e-3, exact_Shat = FALSE,
                                    B_true = NULL) {
  
  N <- nrow(Y); Q <- ncol(Y); P <- ncol(X); J <- max(group)
  stopifnot(all(X[, 1] == 1))
  
  avg_prop <- colMeans(Y / pmax(rowSums(Y), 1))
  mu <- log(avg_prop + 1e-8); mu <- mu - mean(mu)
  phi <- matrix(0, P, Q)
  
  gm <- matrix(0, J, Q)
  for (j in 1:J) {
    idx <- which(group == j)
    if (length(idx) > 0) {
      gs <- colSums(Y[idx, , drop = FALSE])
      gm[j, ] <- log((gs + 1e-5) / (sum(gs) + Q * 1e-5))
    }
  }
  gm_c <- sweep(gm, 2, colMeans(gm))
  sv0 <- svd(t(gm_c), nu = K, nv = K)
  B <- sv0$u %*% diag(pmax(sv0$d[1:K] * 0.5, 0.1), K)
  B <- apply_PLT(B)
  sigma2 <- sigma2_init
  
  log_ev <- numeric(max_iter)
  sigma2_trace <- numeric(max_iter)
  B_dist_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
  converged <- FALSE
  em <- 0
  
  for (em in 1:max_iter) {
    
    es <- estep_wb(J, group, Y, X, M, mu, phi, B, sigma2,
                   max_iter = estep_max_iter, gtol = estep_gtol, exact_Shat = exact_Shat)
    lambda_hat <- es$lambda_hat
    S_hat <- es$S_hat
    
    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma +
      0.5 * es$ld_S_total
    
    mp <- mstep_phi_wb(Y, X, group, lambda_hat, mu, phi,
                       lambda_phi = lambda_phi)
    mu <- mp$mu; phi <- mp$phi
    
    S_obs <- tcrossprod(lambda_hat) / J
    for (j in 1:J) S_obs <- S_obs + S_hat[[j]] / J
    rt <- rubin_thayer_wb(S_obs, K, B_init = B, sigma2_init = sigma2)
    B <- apply_PLT(rt$B)
    sigma2 <- rt$sigma2
    
    sigma2_trace[em] <- sigma2
    if (!is.null(B_true)) B_dist_trace[em] <- subspace_dist(B, B_true)
    
    if (verbose) {
      bd_str <- if (!is.null(B_true)) sprintf("  B_dist = %.4f", B_dist_trace[em]) else ""
      cat(sprintf("iter %3d  log_ev = %.4f  sigma2 = %.4f%s\n",
                  em, log_ev[em], sigma2, bd_str))
    }
    
    if (em > 1 &&
        abs(log_ev[em] - log_ev[em - 1]) <
        tol * (abs(log_ev[em - 1]) + 1)) { converged <- TRUE; break }
  }
  
  list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
       lambda_hat = lambda_hat, S_hat = S_hat,
       log_evidence = log_ev[1:em], converged = converged, iterations = em,
       sigma2_trace = sigma2_trace[1:em],
       B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em] else NULL)
}
