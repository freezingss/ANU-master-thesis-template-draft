# G-structured code with robustness check

# Calculation of likelihood and curvature
group_g_and_W <- function(lambda_j, Y_j, X_j, M_j, mu, phi) {
  Q <- ncol(Y_j)
  fixed <- matrix(mu, nrow(Y_j), Q, byrow = TRUE) + X_j %*% phi
  eta <- sweep(fixed, 2, lambda_j, "+")
  pi_mat <- row_softmax(eta)
  Mpi <- sweep(pi_mat, 1, M_j, "*")
  g_val <- sum(Y_j * eta) - sum(M_j * row_logsumexp(eta))
  W <- diag(as.numeric(colSums(Mpi)), Q) - crossprod(pi_mat, Mpi)   
  list(g = g_val, W = W)
}

# Hessian and invers
exact_H_and_S <- function(lambda_j, Y_j, X_j, M_j, mu, phi, Sinv_mat) {
  gw <- group_g_and_W(lambda_j, Y_j, X_j, M_j, mu, phi)
  Q <- nrow(Sinv_mat)
  H <- Sinv_mat + gw$W
  cH <- tryCatch(chol(H + 1e-10 * diag(Q)), error = function(e) NULL)
  if (is.null(cH)) return(NULL)
  S <- chol2inv(cH)
  logdetH <- 2 * sum(log(diag(cH)))
  list(g = gw$g, W = gw$W, H = H, S = S, logdetH = logdetH)
}

# Core of lambda-block 
Phi_eval <- function(lambda_j, Y_j, X_j, M_j, mu, phi, Sinv_mat) {
  eH <- exact_H_and_S(lambda_j, Y_j, X_j, M_j, mu, phi, Sinv_mat)
  if (is.null(eH)) return(list(Phi = -Inf, S = NULL, logdetH = NA_real_))
  h_val <- eH$g - 0.5 * sum(lambda_j * as.numeric(Sinv_mat %*% lambda_j))
  Phi <- h_val - 0.5 * eH$logdetH
  list(Phi = Phi, S = eH$S, logdetH = eH$logdetH)
}

# Lambda-block
certified_lambda_block <- function(lambda_prev, Y, X, group, M, mu, phi, B, sigma2,
                                    corr_max_rel = 0.5,
                                    estep_max_iter = 100, estep_gtol = 1e-3) {
  Q <- length(mu)
  J <- max(group)
  Sigma <- B %*% t(B) + sigma2 * diag(Q)
  Sinv_mat <- solve(Sigma)
  # logdetSigma <- log(det(Sigma))
  logdetSigma <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  
  # Woodbury approximation does not satisfy condition of log-det, exact_Shat = TRUE necessarily
  es <- estep_wb(J, group, Y, X, M, mu, phi, B, sigma2,
                 max_iter = estep_max_iter, gtol = estep_gtol, exact_Shat = TRUE)
  lambda_hat <- es$lambda_hat
  S_hat_hmode <- es$S_hat
  log_ev <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total  
  
  # mu_{1,j} calculation
  cc <- correct_lambda_edgeworth(lambda_hat, S_hat_hmode, Y, X, group, M, mu, phi)

  lambda_new <- lambda_prev
  S_hat_new <- vector("list", J)
  Phi_total <- 0
  tier_used <- integer(J)   

  for (j in seq_len(J)) {
    idx <- which(group == j)
    
    # Degenerate case (unobserved group)
    if (!length(idx)) { 
      lambda_new[, j] <- lambda_prev[, j]
      S_hat_new[[j]] <- Sigma
      logdetSinv <- -logdetSigma
      Phi_total <- Phi_total - 0.5 * sum(lambda_prev[, j] * as.numeric(Sinv_mat %*% lambda_prev[, j])) -
        0.5 * logdetSinv
      tier_used[j] <- 5L
      next
    }

    Yj <- Y[idx, , drop = FALSE]
    Xj <- X[idx, , drop = FALSE]
    Mj <- M[idx]

    mu1 <- cc$lambda_corrected[, j] - lambda_hat[, j]
    # Threshold
    rel <- cc$rel_size[j]
    if (is.finite(rel) && rel > corr_max_rel) mu1 <- mu1 * (corr_max_rel / rel)

    # Unchanged lambda baseline
    baseline <- Phi_eval(lambda_prev[, j], Yj, Xj, Mj, mu, phi, Sinv_mat)$Phi

    candidates <- list(
      list(lam = lambda_hat[, j] + mu1, tier = 1L),
      list(lam = lambda_hat[, j] + 0.5 * mu1, tier = 2L),
      list(lam = lambda_hat[, j] + 0.25 * mu1, tier = 3L),
      list(lam = lambda_hat[, j], tier = 4L),
      list(lam = lambda_prev[, j], tier = 5L)  
    )

    accepted <- NULL
    for (cand in candidates) {
      ev <- Phi_eval(cand$lam, Yj, Xj, Mj, mu, phi, Sinv_mat)
      if (is.finite(ev$Phi) && ev$Phi >= baseline - 1e-8) {
        accepted <- list(lam = cand$lam, Phi = ev$Phi, S = ev$S, tier = cand$tier)
        break
      }
    }
    if (is.null(accepted)) {
      ev <- Phi_eval(lambda_prev[, j], Yj, Xj, Mj, mu, phi, Sinv_mat)
      accepted <- list(lam = lambda_prev[, j], Phi = ev$Phi, S = ev$S, tier = 5L)
    }

    lambda_new[, j] <- accepted$lam
    S_hat_new[[j]] <- accepted$S
    Phi_total <- Phi_total + accepted$Phi
    tier_used[j] <- accepted$tier
  }

  G_val <- Phi_total - 0.5 * J * logdetSigma

  list(lambda = lambda_new, S_hat = S_hat_new, G = G_val, log_ev = log_ev,
       tier_used = tier_used,
       corr_rel_size = mean(pmin(cc$rel_size, corr_max_rel), na.rm = TRUE))
}

# theta-block
theta_block_update <- function(lambda_cur, S_hat_cur, Y, X, group, M, mu, phi, B, sigma2, K,
                                fix_sigma2 = NULL) {
  Q <- length(mu) 
  J <- max(group)

  S_lam <- tcrossprod(lambda_cur) / J
  S_S <- matrix(0, Q, Q)
  for (j in seq_len(J)) S_S <- S_S + S_hat_cur[[j]] / J
  S_obs_tilde <- S_lam + S_S

  if (is.numeric(fix_sigma2)) {
    sigma2_new <- fix_sigma2
    Ssym <- (S_obs_tilde + t(S_obs_tilde)) / 2
    eS <- eigen(Ssym, symmetric = TRUE)
    lamK <- eS$values[1:K]; Uk <- eS$vectors[, 1:K, drop = FALSE]
    B_new <- apply_PLT(Uk %*% diag(sqrt(pmax(lamK - sigma2_new, 0)), K))
  } else {
    rt <- rubin_thayer_wb(S_obs_tilde, K, B_init = B, sigma2_init = sigma2)
    B_new <- apply_PLT(rt$B); sigma2_new <- rt$sigma2
  }

  obj <- function(mu_, phi_) {
    total <- 0
    for (j in seq_len(J)) {
      idx <- which(group == j)
      if (!length(idx)) next
      gw <- group_g_and_W(lambda_cur[, j], Y[idx, , drop = FALSE], X[idx, , drop = FALSE], M[idx], mu_, phi_)
      total <- total + gw$g - 0.5 * sum(S_hat_cur[[j]] * gw$W)   
    }
    total
  }

  mp <- mstep_phi_wb(Y, X, group, lambda_cur, mu, phi)
  obj_old <- obj(mu, phi)

  accepted <- FALSE
  for (step in c(1, 0.5, 0.25, 0.125)) {
    mu_try  <- (1 - step) * mu  + step * mp$mu
    phi_try <- (1 - step) * phi + step * mp$phi
    if (obj(mu_try, phi_try) >= obj_old - 1e-8) {
      mu_new <- mu_try; phi_new <- phi_try; accepted <- TRUE; break
    }
  }
  if (!accepted) { mu_new <- mu; phi_new <- phi }   

  list(mu = mu_new, phi = phi_new, B = B_new, sigma2 = sigma2_new)
}

fit_pfa_woodbury_G <- function(Y, X, group, K, M = rowSums(Y),
                                max_iter = 80, tol = 1e-8,
                                sigma2_init = 0.3, verbose = FALSE,
                                B_true = NULL, fix_sigma2 = NULL,
                                corr_max_rel = 0.5,
                                estep_max_iter = 100, estep_gtol = 1e-3,
                                B_init = NULL, init_perturb_sd = 0, init_seed = NULL) {
  N <- nrow(Y); Q <- ncol(Y); P <- ncol(X); J <- max(group)
  stopifnot(all(X[, 1] == 1))

  avg_prop <- colMeans(Y / pmax(rowSums(Y), 1))
  mu <- log(avg_prop + 1e-8); mu <- mu - mean(mu)
  phi <- matrix(0, P, Q)

  if (!is.null(B_init)) {
    B <- B_init 
  } else {
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
  }
  if (init_perturb_sd > 0) {
    if (!is.null(init_seed)) set.seed(init_seed)
    B <- B + matrix(rnorm(length(B), sd = init_perturb_sd), nrow(B), ncol(B))
  }
  B <- apply_PLT(B)
  sigma2 <- if (is.numeric(fix_sigma2)) fix_sigma2 else sigma2_init

  es0 <- estep_wb(J, group, Y, X, M, mu, phi, B, sigma2,
                   max_iter = estep_max_iter, gtol = estep_gtol, exact_Shat = TRUE)
  lambda_cur <- es0$lambda_hat
  S_hat_cur <- es0$S_hat

  G_trace <- numeric(max_iter); log_ev_trace <- numeric(max_iter)
  sigma2_trace <- numeric(max_iter)
  B_dist_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
  tier_trace <- matrix(NA_integer_, max_iter, J)
  em <- 0; converged <- FALSE

  for (em in 1:max_iter) {
    tb <- theta_block_update(lambda_cur, S_hat_cur, Y, X, group, M, mu, phi, B, sigma2, K,
                              fix_sigma2 = fix_sigma2)
    mu <- tb$mu; phi <- tb$phi; B <- tb$B; sigma2 <- tb$sigma2

    lb <- certified_lambda_block(lambda_cur, Y, X, group, M, mu, phi, B, sigma2,
                                  corr_max_rel = corr_max_rel,
                                  estep_max_iter = estep_max_iter, estep_gtol = estep_gtol)
    lambda_cur <- lb$lambda; S_hat_cur <- lb$S_hat

    G_trace[em] <- lb$G
    log_ev_trace[em] <- lb$log_ev
    sigma2_trace[em] <- sigma2
    tier_trace[em, ] <- lb$tier_used
    if (!is.null(B_true)) B_dist_trace[em] <- subspace_dist(B, B_true)

    if (em > 1 && G_trace[em] < G_trace[em - 1] - 1e-6) {
      warning(sprintf(
        "iter %d: G decreased (%.6f -> %.6f)",
        em, G_trace[em - 1], G_trace[em]))
    }

    if (verbose) {
      bd <- if (!is.null(B_true)) sprintf("  B_dist=%.4f", B_dist_trace[em]) else ""
      cat(sprintf("iter %3d  G=%.4f  log_ev=%.4f  sigma2=%.4f  corr_sz=%.4f%s\n",
                   em, G_trace[em], log_ev_trace[em], sigma2, lb$corr_rel_size, bd))
    }

    if (em > 1 && abs(G_trace[em] - G_trace[em - 1]) < tol * (abs(G_trace[em - 1]) + 1)) {
      converged <- TRUE; break
    }
  }

  list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
       lambda_hat = lambda_cur, S_hat = S_hat_cur,
       G = G_trace[1:em], log_evidence = log_ev_trace[1:em],
       sigma2_trace = sigma2_trace[1:em],
       B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em] else NULL,
       tier_used = tier_trace[1:em, , drop = FALSE],
       converged = converged, iterations = em)
}

# Fully random multistart for the robustness check

fit_pfa_woodbury_G_multistart <- function(Y, X, group, K, M = rowSums(Y),
                                           n_starts = 5, init_perturb_sd = 0.3,
                                           fully_random_start = FALSE,
                                           base_seed = 1,
                                           max_iter = 80, tol = 1e-8,
                                           sigma2_init = 0.3, verbose = FALSE,
                                           B_true = NULL, fix_sigma2 = NULL,
                                           corr_max_rel = 0.5,
                                           estep_max_iter = 100, estep_gtol = 1e-3) {
  Q <- ncol(Y); J <- max(group)

  gm <- matrix(0, J, Q)
  for (j in 1:J) {
    idx <- which(group == j)
    if (length(idx) > 0) {
      gs <- colSums(Y[idx, , drop = FALSE])
      gm[j, ] <- log((gs + 1e-5) / (sum(gs) + Q * 1e-5))
    }
  }
  gm_c <- sweep(gm, 2, colMeans(gm))
  sv0 <- svd(t(gm_c), nu = 0, nv = 0)
  magnitude <- pmax(sv0$d[1:K] * 0.5, 0.1)

  random_B_init <- function(seed) {
    set.seed(seed)
    Z <- matrix(rnorm(Q * K), Q, K)
    Uk <- qr.Q(qr(Z))[, 1:K, drop = FALSE]   
    Uk %*% diag(magnitude, K)
  }

  fits <- vector("list", n_starts)
  G_final <- numeric(n_starts)
  Bd_final <- if (!is.null(B_true)) numeric(n_starts) else NULL
  iters <- integer(n_starts)
  n_G_decrease_warnings <- integer(n_starts)
  mode_used <- character(n_starts)

  for (s in seq_len(n_starts)) {
    seed_s <- base_seed + s

    if (s == 1) {
      B_init_s <- NULL; perturb_sd_s <- 0; mode_used[s] <- "default"
    } else if (fully_random_start) {
      B_init_s <- random_B_init(seed_s); perturb_sd_s <- 0; mode_used[s] <- "fully_random"
    } else {
      B_init_s <- NULL; perturb_sd_s <- init_perturb_sd; mode_used[s] <- "perturbed"
    }

    warn_count <- 0L
    fit_s <- withCallingHandlers(
      fit_pfa_woodbury_G(Y, X, group, K, M = M,
                          max_iter = max_iter, tol = tol, sigma2_init = sigma2_init,
                          verbose = verbose, B_true = B_true, fix_sigma2 = fix_sigma2,
                          corr_max_rel = corr_max_rel,
                          estep_max_iter = estep_max_iter, estep_gtol = estep_gtol,
                          B_init = B_init_s, init_perturb_sd = perturb_sd_s, init_seed = seed_s),
      warning = function(w) {
        if (grepl("^iter .*: G decreased", conditionMessage(w))) warn_count <<- warn_count + 1L
        invokeRestart("muffleWarning")
      })

    fits[[s]] <- fit_s
    G_final[s] <- tail(fit_s$G, 1)
    iters[s] <- fit_s$iterations
    n_G_decrease_warnings[s] <- warn_count
    if (!is.null(B_true)) Bd_final[s] <- tail(fit_s$B_dist_trace, 1)

    cat(sprintf("start %2d (%-12s seed=%d): iters=%3d  G_final=%.4f%s%s\n",
                s, mode_used[s], seed_s, iters[s], G_final[s],
                if (!is.null(B_true)) sprintf("  Bd_final=%.4f", Bd_final[s]) else "",
                if (warn_count > 0) sprintf("  [%d G-decrease warning(s) -- check for a bug]", warn_count) else ""))
  }

  summary_df <- data.frame(start = seq_len(n_starts), mode = mode_used, iters = iters, G_final = G_final)
  if (!is.null(B_true)) summary_df$Bd_final <- Bd_final
  summary_df <- summary_df[order(-summary_df$G_final), ]

  cat("\n ranking by final G (highest first) \n")
  print(summary_df, row.names = FALSE)
  if (!is.null(B_true)) {
    cat("\n ranking by final B_dist (lowest = best first) \n")
    print(summary_df[order(summary_df$Bd_final), ], row.names = FALSE)
  }

  list(fits = fits, summary = summary_df, best_by_G = fits[[which.max(G_final)]],
       best_by_Bd = if (!is.null(B_true)) fits[[which.min(Bd_final)]] else NULL)
}
