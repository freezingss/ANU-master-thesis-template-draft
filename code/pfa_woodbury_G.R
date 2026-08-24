# =============================================================================
# pfa_woodbury_G.R
#
# Pure function library for the monotone MM algorithm derived in
# mm_monotonicity_fa_dmr.pdf (Section 4: joint objective G, Theorem 4.5 ascent
# guarantee; Section 5: the lambda-block reduces to a CERTIFIED Edgeworth
# step; Section 7: code mapping). No plotting, no comparison/run script --
# call fit_pfa_woodbury_G() from your own scripts.
#
# G(theta, Lambda) = sum_j Phi_j(theta, lambda_j) - (J/2) log|Sigma|,
# Phi_j(theta, lambda) = g_j(lambda; beta) - 0.5 lambda' Sigma^{-1} lambda
#                          - 0.5 log|Sigma^{-1} + W_j(lambda; beta)|
#
# One iteration = (A) theta-block, ascends the EVERYWHERE-TANGENT surrogate Q
#                     built from the log-det tangent-plane inequality
#                     (Lemma 2.1 / Prop 4.4), using the Lambda accepted at the
#                     END of the previous iteration's lambda-block;
#                 (B) lambda-block, re-solves each group's h-mode under the
#                     NEW theta, proposes the Edgeworth-corrected point
#                     (Prop 5.1: this is exactly the Newton step for Phi_j),
#                     and CERTIFIES it -- accepts the first candidate in
#                     {full corr, half, quarter, h-mode, unchanged} that does
#                     not decrease Phi_j relative to the group's old lambda.
#                     "unchanged" always certifies (Phi_j(theta,lambda_old) is
#                     trivially >= itself), so this can never get stuck.
#
# By Theorem 4.5, G is non-decreasing across iterations PROVIDED each block's
# certification is honored. fit_pfa_woodbury_G() asserts this at every
# iteration and warns loudly if it is ever violated -- that would mean a bug,
# not a property of the theory.
#
# sigma2 handling: pass fix_sigma2 = NULL for the unfrozen algorithm, or
# fix_sigma2 = <numeric value> (e.g. a REML-corrected value) to run the
# frozen-sigma2 variant. Both are the SAME function; there is no separate
# "G-frozen" function.
#
# Depends on (source BEFORE this file): pfa_woodbury.R (row_softmax,
# row_logsumexp, laplace_lambda_j_wb2, estep_wb, mstep_phi_wb, rubin_thayer_wb,
# apply_PLT), pfa_woodbury_lambda_corrected.R (correct_lambda_edgeworth),
# basic_functions.R (subspace_dist, apply_PLT if it lives there instead).
#
# NOTE ON COST: unlike the old algorithm, this one needs the EXACT per-group
# curvature W_j = M(diag(pi)-pi pi') (not the diagonal Woodbury majorizer used
# for the old algorithm's Newton DIRECTION), because the log-det tangent-plane
# argument in Prop 4.4 is only valid for the exact Hessian. For Q in the
# range used so far (Q ~ 50) this is a cheap Q x Q Cholesky per evaluation and
# is not a bottleneck; if Q grows much larger, exact_H_and_S() below is the
# one place that would need a Woodbury shortcut re-introduced.
# =============================================================================


# ---------------------------------------------------------------------------
# Exact per-group curvature. Distinct from the Woodbury-approximate curvature
# used internally by laplace_lambda_j_wb2()'s Newton DIRECTION -- that
# approximation is fine for finding a search direction, but Prop 4.4's tangent
# plane needs the exact H_j = Sigma^{-1} + W_j at the point being evaluated.
# ---------------------------------------------------------------------------
group_g_and_W <- function(lambda_j, Y_j, X_j, M_j, mu, phi) {
  Q <- ncol(Y_j)
  fixed <- matrix(mu, nrow(Y_j), Q, byrow = TRUE) + X_j %*% phi
  eta <- sweep(fixed, 2, lambda_j, "+")
  pi_mat <- row_softmax(eta)
  Mpi <- sweep(pi_mat, 1, M_j, "*")
  g_val <- sum(Y_j * eta) - sum(M_j * row_logsumexp(eta))
  W <- diag(as.numeric(colSums(Mpi)), Q) - crossprod(pi_mat, Mpi)   # exact W_j, eq. (2.3) in the note
  list(g = g_val, W = W)
}

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

# Phi_j(theta, lambda) of eq. (3.1) in the note, at an arbitrary lambda
# (not necessarily the h-mode). This is the quantity every certification
# step in Section 5 / the lambda-block compares against.
Phi_eval <- function(lambda_j, Y_j, X_j, M_j, mu, phi, Sinv_mat) {
  eH <- exact_H_and_S(lambda_j, Y_j, X_j, M_j, mu, phi, Sinv_mat)
  if (is.null(eH)) return(list(Phi = -Inf, S = NULL, logdetH = NA_real_))
  h_val <- eH$g - 0.5 * sum(lambda_j * as.numeric(Sinv_mat %*% lambda_j))
  Phi <- h_val - 0.5 * eH$logdetH
  list(Phi = Phi, S = eH$S, logdetH = eH$logdetH)
}


# ---------------------------------------------------------------------------
# (B) The lambda-block: re-solve h-modes under the NEW theta, propose the
# Edgeworth step (Prop 5.1), certify against the OLD lambda under the NEW
# theta (this is exactly what Theorem 4.5's condition (B) requires), and fall
# back gracefully. Also returns the standard log_ev at the fresh h-mode, for
# side-by-side comparison purposes only -- it is not the monitored quantity.
# ---------------------------------------------------------------------------
certified_lambda_block <- function(lambda_prev, Y, X, group, M, mu, phi, B, sigma2,
                                    corr_max_rel = 0.5,
                                    estep_max_iter = 100, estep_gtol = 1e-3) {
  Q <- length(mu); J <- max(group)
  Sigma <- B %*% t(B) + sigma2 * diag(Q)
  Sinv_mat <- solve(Sigma)
  logdetSigma <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)

  # fresh h-modes at the new theta, with EXACT curvature (exact_Shat = TRUE) --
  # needed because the certification below and the reported G both require
  # the exact H_j, not the diagonal Woodbury majorizer.
  es <- estep_wb(J, group, Y, X, M, mu, phi, B, sigma2,
                 max_iter = estep_max_iter, gtol = estep_gtol, exact_Shat = TRUE)
  lambda_hat <- es$lambda_hat
  S_hat_hmode <- es$S_hat
  log_ev <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total  # comparison only, not the monitor

  cc <- correct_lambda_edgeworth(lambda_hat, S_hat_hmode, Y, X, group, M, mu, phi)

  lambda_new <- lambda_prev
  S_hat_new  <- vector("list", J)
  Phi_total  <- 0
  tier_used  <- integer(J)   # 1=full corr, 2=half, 3=quarter, 4=h-mode, 5=unchanged (fallback)

  for (j in seq_len(J)) {
    idx <- which(group == j)

    if (!length(idx)) {   # empty group: no data term, Phi_j has closed form
      lambda_new[, j] <- lambda_prev[, j]
      S_hat_new[[j]]  <- Sigma
      logdetSinv <- -logdetSigma
      Phi_total <- Phi_total - 0.5 * sum(lambda_prev[, j] * as.numeric(Sinv_mat %*% lambda_prev[, j])) -
        0.5 * logdetSinv
      tier_used[j] <- 5L
      next
    }

    Yj <- Y[idx, , drop = FALSE]; Xj <- X[idx, , drop = FALSE]; Mj <- M[idx]

    mu1 <- cc$lambda_corrected[, j] - lambda_hat[, j]
    rel <- cc$rel_size[j]
    if (is.finite(rel) && rel > corr_max_rel) mu1 <- mu1 * (corr_max_rel / rel)

    baseline <- Phi_eval(lambda_prev[, j], Yj, Xj, Mj, mu, phi, Sinv_mat)$Phi

    candidates <- list(
      list(lam = lambda_hat[, j] + mu1,        tier = 1L),
      list(lam = lambda_hat[, j] + 0.5 * mu1,  tier = 2L),
      list(lam = lambda_hat[, j] + 0.25 * mu1, tier = 3L),
      list(lam = lambda_hat[, j],              tier = 4L),
      list(lam = lambda_prev[, j],             tier = 5L)   # always certifies: ties `baseline` exactly
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
    S_hat_new[[j]]  <- accepted$S
    Phi_total <- Phi_total + accepted$Phi
    tier_used[j] <- accepted$tier
  }

  G_val <- Phi_total - 0.5 * J * logdetSigma

  list(lambda = lambda_new, S_hat = S_hat_new, G = G_val, log_ev = log_ev,
       tier_used = tier_used,
       corr_rel_size = mean(pmin(cc$rel_size, corr_max_rel), na.rm = TRUE))
}


# ---------------------------------------------------------------------------
# (A) The theta-block: Sigma-part is the EXACT maximizer of its share of the
# surrogate (Section 4.3, "Sigma-part") -- identical math to the existing
# rubin_thayer_wb / frozen-sigma2 eigen update, only the input S_obs_tilde is
# now built from the ACCEPTED lambda/S_hat pair, not a freshly recomputed
# h-mode. Pass fix_sigma2 = NULL for the unfrozen update, or a numeric value
# to run the frozen-sigma2 variant (same function, same guarantee -- a
# constrained theta-block is still valid per Theorem 4.5 as long as the
# constraint set contains the current iterate). beta-part adds the
# curvature-feedback term of eq. (4.7) that the old mstep_phi_wb ignores; it
# is certified by backtracking toward the old beta, with "no change" as the
# always-valid fallback.
# ---------------------------------------------------------------------------
theta_block_update <- function(lambda_cur, S_hat_cur, Y, X, group, M, mu, phi, B, sigma2, K,
                                fix_sigma2 = NULL) {
  Q <- length(mu); J <- max(group)

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
      total <- total + gw$g - 0.5 * sum(S_hat_cur[[j]] * gw$W)   # eq. (4.7): tr(S W) = sum(S * W), both symmetric
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
  if (!accepted) { mu_new <- mu; phi_new <- phi }   # fallback: "no change" always certifies (A)

  list(mu = mu_new, phi = phi_new, B = B_new, sigma2 = sigma2_new)
}


# ---------------------------------------------------------------------------
# Driver. Monitors G (guaranteed non-decreasing by Theorem 4.5) instead of
# log_ev. Also reports log_ev at each iteration's fresh h-mode purely for
# side-by-side comparison with the old algorithm -- log_ev is NOT the
# quantity the stopping rule or any guarantee applies to here.
#
# fix_sigma2 = NULL      -> unfrozen sigma2 (updated every theta-block)
# fix_sigma2 = <numeric> -> sigma2 frozen at that value for every iteration
# ---------------------------------------------------------------------------
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
    B <- B_init   # caller-supplied starting point (e.g. a multi-start draw, or an oracle start)
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

  # seed Lambda^(0): h-modes + EXACT curvature at the initial theta. This
  # plays the role of Lambda^(t) feeding the FIRST theta-block's surrogate.
  es0 <- estep_wb(J, group, Y, X, M, mu, phi, B, sigma2,
                   max_iter = estep_max_iter, gtol = estep_gtol, exact_Shat = TRUE)
  lambda_cur <- es0$lambda_hat
  S_hat_cur  <- es0$S_hat

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
        "iter %d: G decreased (%.6f -> %.6f). Under Theorem 4.5 this should be impossible -- treat as a bug, not a property of the theory.",
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
