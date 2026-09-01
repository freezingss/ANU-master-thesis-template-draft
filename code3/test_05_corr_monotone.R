source("pfa_woodbury_only.R")
source("pfa_woodbury_lam_corrected.R")

Phi_group <- function(Y_j, X_j, M_j, mu, phi, B, sigma2,
                      M_K, M_K_inv, log_det_MK, lambda, want_S = FALSE) {
  N_j <- nrow(Y_j)
  Q <- ncol(Y_j)
  K <- ncol(B)
  fixed <- matrix(rep(mu, each = N_j), N_j, Q) + X_j %*% phi
  Sinv_lam <- lambda / sigma2 - B %*% (M_K_inv %*% (t(B) %*% lambda)) / sigma2^2
  eta <- sweep(fixed, 2, lambda, "+")
  h_val <- sum(Y_j * eta) - sum(M_j * row_logsumexp(eta)) - 0.5 * sum(lambda * Sinv_lam)
  pi_mat <- row_softmax(eta)
  Mpi <- sweep(pi_mat, 1, M_j, "*")
  d <- as.numeric(colSums(Mpi))
  W <- sqrt(M_j) * pi_mat
  dt <- d + 1 / sigma2
  idt <- 1 / dt
  Widt <- sweep(W, 2, idt, "*")
  U <- cbind(t(W), B)
  Delta_11 <- diag(N_j) - Widt %*% t(W)
  Delta_12 <- -Widt %*% B
  Delta_22 <- sigma2^2 * M_K - crossprod(B, idt * B)
  Delta <- rbind(cbind(Delta_11, Delta_12), cbind(t(Delta_12), Delta_22))
  L <- tryCatch(chol(Delta + 1e-10 * diag(nrow(Delta))), error = function(e) NULL)
  if (is.null(L)) return(NULL)
  log_det_Delta <- 2 * sum(log(diag(L)))
  log_det_shat <- -sum(log(dt)) + log_det_MK + 2 * K * log(sigma2) - log_det_Delta
  out <- list(Phi = h_val + 0.5 * log_det_shat, h = h_val, log_det_shat = log_det_shat)
  if (want_S) {
    DeltaInv <- chol2inv(L)
    UdiagIdt <- idt * U
    out$S <- diag(idt, Q) + UdiagIdt %*% DeltaInv %*% t(UdiagIdt)
  }
  out
}

G_eval <- function(J, group, Y, X, M, mu, phi, B, sigma2, Lambda) {
  Q <- length(mu)
  K <- ncol(B)
  M_K <- diag(K) + crossprod(B) / sigma2
  M_K_inv <- tryCatch(solve(M_K), error = function(e) NULL)
  if (is.null(M_K_inv)) return(NULL)
  log_det_MK <- 2 * sum(log(diag(chol(M_K))))
  ldSigma <- Q * log(sigma2) + log_det_MK
  G_val <- -0.5 * J * ldSigma
  for (j in 1:J) {
    idx <- which(group == j)
    if (!length(idx)) {
      G_val <- G_val + 0.5 * ldSigma
      next
    }
    pe <- Phi_group(Y[idx, , drop = FALSE], X[idx, , drop = FALSE], M[idx],
                    mu, phi, B, sigma2, M_K, M_K_inv, log_det_MK, Lambda[, j])
    if (is.null(pe)) return(NULL)
    G_val <- G_val + pe$Phi
  }
  G_val
}

referee_theta <- function(J, group, Y, X, M, cur, cand, Lambda, G_cur, bt_max = 6) {
  s <- 1
  for (bt in 1:bt_max) {
    mu_s <- cur$mu + s * (cand$mu - cur$mu)
    phi_s <- cur$phi + s * (cand$phi - cur$phi)
    B_s <- cur$B + s * (cand$B - cur$B)
    sig_s <- cur$sigma2 + s * (cand$sigma2 - cur$sigma2)
    G_s <- G_eval(J, group, Y, X, M, mu_s, phi_s, B_s, sig_s, Lambda)
    if (!is.null(G_s) && G_s >= G_cur) {
      return(list(mu = mu_s, phi = phi_s, B = B_s, sigma2 = sig_s,
                  G = G_s, step = s, accepted = TRUE))
    }
    s <- s / 2
  }
  list(mu = cur$mu, phi = cur$phi, B = cur$B, sigma2 = cur$sigma2,
       G = G_cur, step = 0, accepted = FALSE)
}

lambda_block <- function(J, group, Y, X, M, mu, phi, B, sigma2,
                         lambda_hat, S_hat, Lambda_prev, want_S = FALSE) {
  Q <- nrow(lambda_hat)
  K <- ncol(B)
  M_K <- diag(K) + crossprod(B) / sigma2
  M_K_inv <- solve(M_K)
  log_det_MK <- 2 * sum(log(diag(chol(M_K))))
  ldSigma <- Q * log(sigma2) + log_det_MK
  cc <- correct_lambda_edgeworth(Q, J, lambda_hat, S_hat, Y, X, group, M, mu, phi)
  mu1 <- cc$lambda_corrected - lambda_hat
  Lambda_new <- lambda_hat
  S_out <- vector("list", J)
  Phi_sum <- 0
  win <- rep(NA_integer_, J)
  for (j in 1:J) {
    idx <- which(group == j)
    if (!length(idx)) {
      Lambda_new[, j] <- 0
      S_out[[j]] <- B %*% t(B) + sigma2 * diag(Q)
      Phi_sum <- Phi_sum + 0.5 * ldSigma
      next
    }
    Yj <- Y[idx, , drop = FALSE]
    Xj <- X[idx, , drop = FALSE]
    Mj <- M[idx]
    cands <- list(lambda_hat[, j] + mu1[, j],
                  lambda_hat[, j] + 0.5 * mu1[, j],
                  lambda_hat[, j],
                  Lambda_prev[, j])
    best_val <- -Inf
    best_id <- 3L
    for (k in seq_along(cands)) {
      pe <- Phi_group(Yj, Xj, Mj, mu, phi, B, sigma2, M_K, M_K_inv, log_det_MK, cands[[k]])
      if (!is.null(pe) && pe$Phi > best_val) {
        best_val <- pe$Phi
        best_id <- k
      }
    }
    if (!is.finite(best_val)) stop("Phi evaluation failed for all candidates in group ", j)
    Lambda_new[, j] <- cands[[best_id]]
    win[j] <- best_id
    Phi_sum <- Phi_sum + best_val
    if (want_S) {
      pw <- Phi_group(Yj, Xj, Mj, mu, phi, B, sigma2, M_K, M_K_inv, log_det_MK,
                      Lambda_new[, j], want_S = TRUE)
      S_out[[j]] <- pw$S
    } else {
      S_out[[j]] <- S_hat[[j]]
    }
  }
  list(Lambda = Lambda_new, S_list = S_out, G = Phi_sum - 0.5 * J * ldSigma,
       win = win, rel_size = cc$rel_size)
}

fit_pfa_G <- function(Y, X, group, K,
                      M = rowSums(Y),
                      max_iter = 200, tol = 1e-4,
                      lambda_phi = 0, sigma2_init = 0.3,
                      verbose = TRUE,
                      estep_max_iter = 100, estep_gtol = 1e-3,
                      B_true = NULL, use_Sbar = FALSE, bt_max = 6) {
  N <- nrow(Y)
  Q <- ncol(Y)
  P <- ncol(X)
  J <- max(group)
  stopifnot(all(X[, 1] == 1))
  avg_prop <- colMeans(Y / pmax(rowSums(Y), 1))
  mu <- log(avg_prop + 1e-8)
  mu <- mu - mean(mu)
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

  es <- estep_wbonly(J, group, Y, X, M, mu, phi, B, sigma2, estep_max_iter, estep_gtol)
  log_ev_cur <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total
  lb <- lambda_block(J, group, Y, X, M, mu, phi, B, sigma2,
                     es$lambda_hat, es$S_hat, es$lambda_hat, want_S = use_Sbar)
  Lambda <- lb$Lambda
  S_list <- if (use_Sbar) lb$S_list else es$S_hat
  G_cur <- lb$G

  G_trace <- numeric(max_iter)
  log_ev_trace <- numeric(max_iter)
  gap_trace <- numeric(max_iter)
  sigma2_trace <- numeric(max_iter)
  B_dist_trace <- if (!is.null(B_true)) numeric(max_iter) else NULL
  step_muphi <- numeric(max_iter)
  step_Bsig <- numeric(max_iter)
  win_counts <- matrix(0L, max_iter, 4)
  mono_violation <- 0L
  converged <- FALSE
  em <- 0

  for (em in 1:max_iter) {
    G_start <- G_cur

    mp <- mstep_phi_wbonly(Y, X, group, Lambda, mu, phi, lambda_phi = lambda_phi)
    r1 <- referee_theta(J, group, Y, X, M,
                        list(mu = mu, phi = phi, B = B, sigma2 = sigma2),
                        list(mu = mp$mu, phi = mp$phi, B = B, sigma2 = sigma2),
                        Lambda, G_cur, bt_max)
    mu <- r1$mu
    phi <- r1$phi
    G_cur <- r1$G
    step_muphi[em] <- r1$step

    S_obs <- tcrossprod(Lambda) / J
    for (j in 1:J) S_obs <- S_obs + S_list[[j]] / J
    rt <- rubin_thayer_wbonly(S_obs, K, B_init = B, sigma2_init = sigma2)
    r2 <- referee_theta(J, group, Y, X, M,
                        list(mu = mu, phi = phi, B = B, sigma2 = sigma2),
                        list(mu = mu, phi = phi, B = apply_PLT(rt$B), sigma2 = rt$sigma2),
                        Lambda, G_cur, bt_max)
    B <- r2$B
    sigma2 <- r2$sigma2
    G_cur <- r2$G
    step_Bsig[em] <- r2$step

    es <- estep_wbonly(J, group, Y, X, M, mu, phi, B, sigma2, estep_max_iter, estep_gtol)
    log_ev_cur <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total

    lb <- lambda_block(J, group, Y, X, M, mu, phi, B, sigma2,
                       es$lambda_hat, es$S_hat, Lambda, want_S = use_Sbar)
    Lambda <- lb$Lambda
    S_list <- if (use_Sbar) lb$S_list else es$S_hat
    G_new <- lb$G
    if (G_new < G_cur - 1e-8) mono_violation <- mono_violation + 1L
    if (G_new < G_start - 1e-8) mono_violation <- mono_violation + 1L

    G_trace[em] <- G_new
    log_ev_trace[em] <- log_ev_cur
    gap_trace[em] <- G_new - log_ev_cur
    sigma2_trace[em] <- sigma2
    if (!is.null(B_true)) B_dist_trace[em] <- subspace_dist(B, B_true)
    win_counts[em, ] <- as.integer(table(factor(lb$win, levels = 1:4)))

    if (verbose) {
      bd <- if (!is.null(B_true)) sprintf("  B_dist=%.4f", B_dist_trace[em]) else ""
      cat(sprintf("iter %3d  G=%.4f  log_ev=%.4f  gap=%.4f  sigma2=%.4f  s1=%.2f  s2=%.2f%s\n",
                  em, G_new, log_ev_cur, gap_trace[em], sigma2,
                  step_muphi[em], step_Bsig[em], bd))
    }

    if (em > 1 && abs(G_new - G_trace[em - 1]) < tol) {
      G_cur <- G_new
      converged <- TRUE
      break
    }
    G_cur <- G_new
  }

  list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
       Lambda = Lambda, lambda_hat = es$lambda_hat, S_hat = es$S_hat,
       G_trace = G_trace[1:em], log_ev_trace = log_ev_trace[1:em],
       gap_trace = gap_trace[1:em], sigma2_trace = sigma2_trace[1:em],
       B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em] else NULL,
       step_muphi = step_muphi[1:em], step_Bsig = step_Bsig[1:em],
       win_counts = win_counts[1:em, , drop = FALSE],
       mono_violation = mono_violation,
       converged = converged, iterations = em)
}

# Small test: J = Q = 50, seed = 1
source("basic_functions.R")
source("sim_data.R")
# source("pfa_G_monotone.R")

check_monotone <- function(x, eps = 1e-6) {
  d <- diff(x)
  viol <- which(d < -eps)
  list(n = length(viol), total = length(d),
       worst = if (length(viol)) min(d) else 0,
       at = viol)
}

J <- 100
Q <- 150
K <- 2
seed <- 1
sigma2_true <- 0.3

dat <- simulate_pfa_data(Q = Q, K = K, J = J, N_per_group = 15,
                         P = 3, sigma2 = sigma2_true, M_rate = 150,
                         seed = seed)

t_corr <- system.time(
  fit_corr <- fit_pfa_woodbury_lam_corr(dat$Y, dat$X, dat$group, K,
                                        max_iter = 80, tol = 1e-4,
                                        verbose = TRUE,
                                        B_true = dat$true$B,
                                        use_lambda_correction = TRUE)
)[1]

t_G <- system.time(
  fit_G <- fit_pfa_G(dat$Y, dat$X, dat$group, K,
                     max_iter = 80, tol = 1e-4,
                     verbose = TRUE,
                     B_true = dat$true$B)
)[1]

mc_corr <- check_monotone(fit_corr$log_evidence)
mc_G_ev <- check_monotone(fit_G$log_ev_trace)
mc_G <- check_monotone(fit_G$G_trace)

cat("\n monotonicity \n")
cat(sprintf("corr  log_ev : %d / %d violations, worst drop %.6f\n",
            mc_corr$n, mc_corr$total, mc_corr$worst))
cat(sprintf("G     log_ev : %d / %d violations, worst drop %.6f\n",
            mc_G_ev$n, mc_G_ev$total, mc_G_ev$worst))
cat(sprintf("G     G      : %d / %d violations, worst drop %.6f\n",
            mc_G$n, mc_G$total, mc_G$worst))
cat(sprintf("G internal mono_violation counter: %d\n", fit_G$mono_violation))
cat(sprintf("gap = G - log_ev, min over run: %.6f (must be >= 0)\n",
            min(fit_G$gap_trace)))

cat("\n estimates \n")
cat(sprintf("corr : iters=%3d  converged=%s  B_dist=%.4f  sigma2=%.4f  time=%.1fs\n",
            fit_corr$iterations, fit_corr$converged,
            tail(fit_corr$B_dist_trace, 1), fit_corr$sigma2, t_corr))
cat(sprintf("G    : iters=%3d  converged=%s  B_dist=%.4f  sigma2=%.4f  time=%.1fs\n",
            fit_G$iterations, fit_G$converged,
            tail(fit_G$B_dist_trace, 1), fit_G$sigma2, t_G))
cat(sprintf("sigma2 true = %.4f\n", sigma2_true))

cat("\n G diagnostics \n")
cat("lambda-block candidate wins (full mu1 / half mu1 / mode / previous):\n")
print(colSums(fit_G$win_counts))
cat(sprintf("mu-phi block: mean step %.3f, damped in %d / %d iters\n",
            mean(fit_G$step_muphi), sum(fit_G$step_muphi < 1), fit_G$iterations))
cat(sprintf("B-sigma block: mean step %.3f, damped in %d / %d iters\n",
            mean(fit_G$step_Bsig), sum(fit_G$step_Bsig < 1), fit_G$iterations))

pdf(sprintf("test_G_monotone_J%d_Q%d_seed%d.pdf", J, Q, seed), width = 11, height = 4.5)
par(mfrow = c(1, 3))
d1 <- c(NA, diff(fit_corr$log_evidence))
plot(fit_corr$log_evidence, type = "l", col = "steelblue",
     xlab = "iteration", ylab = "log_ev",
     main = sprintf("corr log_ev (%d violations)", mc_corr$n))
points(which(d1 < -1e-6), fit_corr$log_evidence[which(d1 < -1e-6)],
       col = "red", pch = 19)
d2 <- c(NA, diff(fit_G$log_ev_trace))
plot(fit_G$log_ev_trace, type = "l", col = "steelblue",
     xlab = "iteration", ylab = "log_ev",
     main = sprintf("G log_ev (%d violations)", mc_G_ev$n))
points(which(d2 < -1e-6), fit_G$log_ev_trace[which(d2 < -1e-6)],
       col = "red", pch = 19)
d3 <- c(NA, diff(fit_G$G_trace))
plot(fit_G$G_trace, type = "l", col = "darkgreen",
     xlab = "iteration", ylab = "G",
     main = sprintf("G certified trace (%d violations)", mc_G$n))
points(which(d3 < -1e-6), fit_G$G_trace[which(d3 < -1e-6)],
       col = "red", pch = 19)
dev.off()

fit_pfa_woodbury_lam_corr_decoupled <- function(Y, X, group, K,
                                                M = rowSums(Y),
                                                max_iter = 60, tol = 1e-4,
                                                lambda_phi = 0, sigma2_init = 0.3,
                                                verbose = FALSE,
                                                estep_max_iter = 100, estep_gtol = 1e-3,
                                                B_true = NULL,
                                                use_lambda_correction = FALSE,
                                                decouple_muphi_lambda = TRUE,
                                                corr_max_rel = 0.5,
                                                fix_sigma2 = NULL,
                                                freeze_tol = 1e-5, freeze_patience = 3,
                                                use_aitken = FALSE, aitken_window = 3, aitken_tol = 1e-4) {
  
  N <- nrow(Y)
  Q <- ncol(Y)
  P <- ncol(X)
  J <- max(group)
  stopifnot(all(X[, 1] == 1))
  stopifnot(is.null(fix_sigma2) || identical(fix_sigma2, "auto") || is.numeric(fix_sigma2))
  
  avg_prop <- colMeans(Y / pmax(rowSums(Y), 1))
  mu <- log(avg_prop + 1e-8)
  mu <- mu - mean(mu)
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
  corr_rel_size_trace <- numeric(max_iter)
  corr_n_capped_trace <- integer(max_iter)
  Shat_share_trace <- numeric(max_iter)
  converged <- FALSE
  em <- 0
  
  sigma2_is_frozen <- FALSE
  sigma2_freeze_iter <- NA_integer_
  sigma2_freeze_value <- NA_real_
  freeze_stable_count <- 0
  if (is.numeric(fix_sigma2)) {
    sigma2 <- fix_sigma2
    sigma2_is_frozen <- TRUE
    sigma2_freeze_iter <- 0L
    sigma2_freeze_value <- fix_sigma2
    if (verbose) message(sprintf("sigma2 frozen from iter 0 at %.4f", fix_sigma2))
  }
  aitken_extrapolate <- function(x3) {
    d <- x3[3] - 2 * x3[2] + x3[1]
    if (!is.finite(d) || abs(d) < 1e-12) return(NA_real_)
    x3[3] - (x3[3] - x3[2])^2 / d
  }
  
  lambda_mode <- NULL
  lambda_corr <- NULL
  
  for (em in 1:max_iter) {
    
    es <- estep_wbonly(J, group, Y, X, M, mu, phi, B, sigma2,
                       estep_max_iter, estep_gtol)
    lambda_mode <- es$lambda_hat
    S_hat <- es$S_hat
    
    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma + 0.5 * es$ld_S_total
    
    n_capped <- 0L
    if (use_lambda_correction) {
      cc <- correct_lambda_edgeworth(Q, J, lambda_mode, S_hat, Y, X, group, M, mu, phi)
      lam_c <- cc$lambda_corrected
      for (j in seq_len(J)) {
        if (is.finite(cc$rel_size[j]) && cc$rel_size[j] > corr_max_rel) {
          mu1 <- lam_c[, j] - lambda_mode[, j]
          lam_c[, j] <- lambda_mode[, j] + mu1 * (corr_max_rel / cc$rel_size[j])
          n_capped <- n_capped + 1L
        }
      }
      lambda_corr <- lam_c
      corr_rel_size_trace[em] <- mean(pmin(cc$rel_size, corr_max_rel), na.rm = TRUE)
    } else {
      lambda_corr <- lambda_mode
      corr_rel_size_trace[em] <- 0
    }
    corr_n_capped_trace[em] <- n_capped
    
    lambda_for_muphi <- if (decouple_muphi_lambda) lambda_mode else lambda_corr
    lambda_for_S <- lambda_corr
    
    mp <- mstep_phi_wbonly(Y, X, group, lambda_for_muphi, mu, phi, lambda_phi = lambda_phi)
    mu <- mp$mu
    phi <- mp$phi
    
    S_lam <- tcrossprod(lambda_for_S) / J
    S_S <- matrix(0, Q, Q)
    for (j in 1:J) S_S <- S_S + S_hat[[j]] / J
    S_obs <- S_lam + S_S
    Shat_share_trace[em] <- sum(diag(S_S)) / sum(diag(S_obs))
    
    if (sigma2_is_frozen) {
      Ssym <- (S_obs + t(S_obs)) / 2
      eS <- eigen(Ssym, symmetric = TRUE)
      lamK <- eS$values[1:K]
      Uk <- eS$vectors[, 1:K, drop = FALSE]
      if (any(lamK < sigma2) && verbose)
        message(sprintf("iter %d: %d factor eigenvalue(s) below frozen sigma2 - clipped", em, sum(lamK < sigma2)))
      B <- apply_PLT(Uk %*% diag(sqrt(pmax(lamK - sigma2, 0)), K))
    } else {
      rt <- rubin_thayer_wbonly(S_obs, K, B_init = B, sigma2_init = sigma2)
      B <- apply_PLT(rt$B)
      sigma2 <- rt$sigma2
    }
    
    sigma2_trace[em] <- sigma2
    if (!is.null(B_true)) B_dist_trace[em] <- subspace_dist(B, B_true)
    
    if (verbose) {
      bd <- if (!is.null(B_true)) sprintf("  B_dist=%.4f", B_dist_trace[em]) else ""
      cs <- if (use_lambda_correction) sprintf("  corr_sz=%.4f capped=%d", corr_rel_size_trace[em], n_capped) else ""
      fz <- if (sigma2_is_frozen) sprintf("  [FROZEN@%.4f]", sigma2_freeze_value) else ""
      cat(sprintf("iter %3d  log_ev=%.4f  sigma2=%.4f  Shat_share=%.3f%s%s%s\n",
                  em, log_ev[em], sigma2, Shat_share_trace[em], bd, cs, fz))
    }
    
    if (identical(fix_sigma2, "auto") && !sigma2_is_frozen) {
      if (use_aitken && em >= aitken_window + 1) {
        est_new <- aitken_extrapolate(sigma2_trace[(em - aitken_window + 1):em])
        est_old <- aitken_extrapolate(sigma2_trace[(em - aitken_window):(em - 1)])
        if (is.finite(est_new) && is.finite(est_old) &&
            abs(est_new - est_old) < aitken_tol * (abs(est_old) + 1e-8)) {
          sigma2 <- est_new
          sigma2_is_frozen <- TRUE; sigma2_freeze_iter <- em; sigma2_freeze_value <- est_new
          if (verbose) message(sprintf("iter %d: Aitken freeze at %.4f", em, est_new))
        }
      } else if (!use_aitken && em >= 2) {
        if (abs(sigma2_trace[em] - sigma2_trace[em - 1]) < freeze_tol) {
          freeze_stable_count <- freeze_stable_count + 1
        } else freeze_stable_count <- 0
        if (freeze_stable_count >= freeze_patience) {
          sigma2_is_frozen <- TRUE; sigma2_freeze_iter <- em; sigma2_freeze_value <- sigma2_trace[em]
          if (verbose) message(sprintf("iter %d: tol freeze at %.4f", em, sigma2_freeze_value))
        }
      }
    }
    
    if (em > 1 &&
        abs(log_ev[em] - log_ev[em - 1]) < tol) {
      converged <- TRUE; break
    }
  }
  
  list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
       lambda_hat = lambda_corr, lambda_mode = lambda_mode, lambda_corrected = lambda_corr,
       S_hat = S_hat,
       log_evidence = log_ev[1:em], converged = converged, iterations = em,
       sigma2_trace = sigma2_trace[1:em],
       B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em] else NULL,
       corr_rel_size_trace = corr_rel_size_trace[1:em],
       corr_n_capped_trace = corr_n_capped_trace[1:em],
       Shat_share_trace = Shat_share_trace[1:em],
       sigma2_freeze_iter = sigma2_freeze_iter,
       sigma2_freeze_value = sigma2_freeze_value,
       decoupled = decouple_muphi_lambda)
}

fit_coupled <- fit_pfa_woodbury_lam_corr(dat$Y, dat$X, dat$group, K,
                                         max_iter = 80, tol = 1e-4,
                                         verbose = FALSE,
                                         B_true = dat$true$B,
                                         use_lambda_correction = TRUE)

fit_decoupled <- fit_pfa_woodbury_lam_corr_decoupled(dat$Y, dat$X, dat$group, K,
                                                     max_iter = 80, tol = 1e-4,
                                                     verbose = FALSE,
                                                     B_true = dat$true$B,
                                                     use_lambda_correction = TRUE,
                                                     decouple_muphi_lambda = TRUE)

mc_coupled <- check_monotone(fit_coupled$log_evidence)
mc_decoupled <- check_monotone(fit_decoupled$log_evidence)

cat(sprintf("coupled    log_ev : %d / %d violations, worst drop %.6f\n",
            mc_coupled$n, mc_coupled$total, mc_coupled$worst))
cat(sprintf("decoupled  log_ev : %d / %d violations, worst drop %.6f\n",
            mc_decoupled$n, mc_decoupled$total, mc_decoupled$worst))
cat(sprintf("coupled   : B_dist=%.4f  sigma2=%.4f  iters=%d  converged=%s\n",
            tail(fit_coupled$B_dist_trace, 1), fit_coupled$sigma2,
            fit_coupled$iterations, fit_coupled$converged))
cat(sprintf("decoupled : B_dist=%.4f  sigma2=%.4f  iters=%d  converged=%s\n",
            tail(fit_decoupled$B_dist_trace, 1), fit_decoupled$sigma2,
            fit_decoupled$iterations, fit_decoupled$converged))

par(mfrow = c(1, 2))
d1 <- c(NA, diff(fit_coupled$log_evidence))
plot(fit_coupled$log_evidence, type = "l", col = "steelblue",
     xlab = "iteration", ylab = "log_ev",
     main = sprintf("coupled (%d violations)", mc_coupled$n))
points(which(d1 < -1e-6), fit_coupled$log_evidence[which(d1 < -1e-6)], col = "red", pch = 19)
d2 <- c(NA, diff(fit_decoupled$log_evidence))
plot(fit_decoupled$log_evidence, type = "l", col = "steelblue",
     xlab = "iteration", ylab = "log_ev",
     main = sprintf("decoupled (%d violations)", mc_decoupled$n))
points(which(d2 < -1e-6), fit_decoupled$log_evidence[which(d2 < -1e-6)], col = "red", pch = 19)
