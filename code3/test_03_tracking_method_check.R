# Compare 3 methods: base (multinomial) & Woodbury & Woodbury + arjimo backtracking

library(bench)
library(tidyverse)

source("pfa_base.R")
source("pfa_woodbury.R")
source("basic_functions.R")
source("sim_data.R")
source("pfa_woodbury_traced.R")

laplace_lambda_j_wb2_fixedstep <- function(Y_j, X_j, M_j, mu, phi, B, sigma2,
                                           BtB, M_K, M_K_inv, log_det_MK,
                                           max_iter = 100, gtol = 1e-3, exact_Shat = FALSE) {
  N_j <- nrow(Y_j)
  Q <- ncol(Y_j)
  K <- ncol(B)
  fixed <- matrix(rep(mu, each = N_j), N_j, Q) + X_j %*% phi
  Sinv <- function(v) v / sigma2 - B %*% (M_K_inv %*% (t(B) %*% v)) / sigma2^2
  lp <- function(a) {
    eta <- sweep(fixed, 2, a, "+")
    sum(Y_j * eta) - sum(M_j * row_logsumexp(eta)) - 0.5 * sum(a * Sinv(a))
  }
  wb_solve <- function(g, d) {
    dt <- d + 1 / sigma2
    idt <- 1 / dt
    Gam <- t(B) %*% (idt * B) - sigma2^2 * M_K
    L <- tryCatch(chol(-Gam + 1e-10 * diag(K)), error = function(e) NULL)
    if (is.null(L)) return(g * 0.01)
    idt * g - idt * (B %*% backsolve(L, forwardsolve(t(L), -crossprod(B, idt * g))))
  }
  lambda <- rep(0, Q)
  lp_cur <- lp(lambda)
  g <- rep(Inf, Q)
  d <- rep(0, Q)
  for (iter in 1:max_iter) {
    eta <- sweep(fixed, 2, lambda, "+")
    pi <- row_softmax(eta)
    Mpi <- sweep(pi, 1, M_j, "*")
    d <- as.numeric(colSums(Mpi))
    g <- as.numeric(colSums(Y_j - Mpi)) - Sinv(lambda)
    if (max(abs(g)) < gtol) break
    dir <- wb_solve(g, d)
    lambda <- lambda + dir
    lp_cur <- lp(lambda)
  }
  eta <- sweep(fixed, 2, lambda, "+")
  pi <- row_softmax(eta)
  Mpi <- sweep(pi, 1, M_j, "*")
  d <- as.numeric(colSums(Mpi))
  g <- as.numeric(colSums(Y_j - Mpi)) - Sinv(lambda)
  if (exact_Shat) {
    Sinv_mat <- diag(1 / sigma2, Q) - B %*% M_K_inv %*% t(B) / sigma2^2
    negH <- diag(as.numeric(colSums(Mpi)), Q) - crossprod(pi, Mpi) + Sinv_mat
    cL <- chol(negH)
    S_hat <- chol2inv(cL)
    log_det_shat <- -2 * sum(log(diag(cL)))
  } else {
    dt  <- d + 1 / sigma2
    idt <- 1 / dt
    Gam <- t(B) %*% (idt * B) - sigma2^2 * M_K
    L <- chol(-Gam + 1e-10 * diag(K))
    log_det_shat <- -(sum(log(dt)) - 2 * K * log(sigma2) - log_det_MK +
                        2 * sum(log(diag(L))))
    GamInv <- chol2inv(L)
    IdtB <- idt * B
    S_hat <- diag(idt, Q) + IdtB %*% GamInv %*% t(IdtB)
  }
  list(lambda_hat = lambda, S_hat = S_hat, lp_mode = lp_cur,
       log_det_shat = log_det_shat, n_iter = iter, grad_norm = max(abs(g)))
}

estep_wb <- function(J, group, Y, X, M, mu, phi, B, sigma2,
                     max_iter = 100, gtol = 1e-3, exact_Shat = FALSE, armijo_bt = TRUE) {
  Q <- length(mu); K <- ncol(B)
  
  BtB <- crossprod(B)
  M_K <- diag(K) + BtB / sigma2
  M_K_inv <- solve(M_K)
  log_det_MK <- 2 * sum(log(diag(chol(M_K))))
  ldSigma <- Q * log(sigma2) + log_det_MK
  
  lambda_hat <- matrix(0, Q, J); S_hat <- vector("list", J)
  lp_total <- 0; ld_S_total <- 0
  n_iter_vec <- numeric(J); grad_norm_vec <- numeric(J)
  
  mode_fun <- if (armijo_bt) laplace_lambda_j_wb2 else laplace_lambda_j_wb2_fixedstep
  
  for (j in 1:J) {
    idx <- which(group == j)
    if (!length(idx)) {
      S_hat[[j]] <- B %*% t(B) + sigma2 * diag(Q)
      ld_S_total <- ld_S_total + ldSigma
      n_iter_vec[j] <- NA; grad_norm_vec[j] <- NA
      next
    }
    res <- mode_fun(Y[idx, , drop = FALSE], X[idx, , drop = FALSE],
                    M[idx], mu, phi, B, sigma2,
                    BtB = BtB, M_K = M_K, M_K_inv = M_K_inv,
                    log_det_MK = log_det_MK,
                    max_iter = max_iter, gtol = gtol, exact_Shat = exact_Shat)
    lambda_hat[, j] <- res$lambda_hat
    S_hat[[j]] <- res$S_hat
    lp_total <- lp_total + res$lp_mode
    ld_S_total <- ld_S_total + res$log_det_shat
    n_iter_vec[j] <- res$n_iter
    grad_norm_vec[j] <- res$grad_norm
  }
  list(lambda_hat = lambda_hat, S_hat = S_hat, lp_total = lp_total,
       ld_S_total = ld_S_total, log_det_Sigma = ldSigma,
       n_iter_vec = n_iter_vec, grad_norm_vec = grad_norm_vec)
}

fit_pfa_woodbury_traced <- function(Y, X, group, K,
                                    M = rowSums(Y),
                                    max_iter = 60, tol = 1e-4,
                                    lambda_phi = 0, sigma2_init = 0.3,
                                    verbose = FALSE,
                                    estep_max_iter = 100, estep_gtol = 1e-3, exact_Shat = FALSE,
                                    armijo_bt = TRUE,
                                    B_true = NULL) {
  
  N <- nrow(Y)
  Q <- ncol(Y)
  P <- ncol(X)
  J <- max(group)
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
  es <- NULL
  
  for (em in 1:max_iter) {
    
    es <- estep_wb(J, group, Y, X, M, mu, phi, B, sigma2,
                   max_iter = estep_max_iter, gtol = estep_gtol, exact_Shat = exact_Shat,
                   armijo_bt = armijo_bt)
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
        abs(log_ev[em] - log_ev[em - 1]) < tol) { converged <- TRUE; break }
  }
  
  list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
       lambda_hat = lambda_hat, S_hat = S_hat,
       log_evidence = log_ev[1:em], converged = converged, iterations = em,
       sigma2_trace = sigma2_trace[1:em],
       B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em] else NULL,
       inner_iter_mean_final = mean(es$n_iter_vec, na.rm = TRUE),
       grad_norm_mean_final = mean(es$grad_norm_vec, na.rm = TRUE),
       inner_maxiter_hit_pct_final = mean(es$n_iter_vec >= estep_max_iter, na.rm = TRUE))
}

estep_dense <- function(J, group, Y, X, M, mu, phi, B, sigma2,
                        max_iter = 100, gtol = 1e-3) {
  Q <- length(mu)
  Sigma <- B %*% t(B) + sigma2 * diag(Q)
  Sigma_inv <- tryCatch({
    cholSigma <- Cholesky(Sigma)
    solve(cholSigma, Diagonal(Q))
  }, error = function(e) {diag(Q) / sigma2})
  ldSigma <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  
  lambda_hat <- matrix(0, Q, J)
  S_hat <- vector("list", J)
  lp_total <- 0
  ld_S_total <- 0
  n_iter_vec <- numeric(J); grad_norm_vec <- numeric(J)
  
  for (j in 1:J) {
    idx <- which(group == j)
    if (!length(idx)) {
      S_hat[[j]] <- Sigma
      ld_S_total <- ld_S_total + ldSigma
      n_iter_vec[j] <- NA; grad_norm_vec[j] <- NA
      next
    }
    res <- laplace_lambda_j_dense(Y[idx, , drop = FALSE], X[idx, , drop = FALSE],
                                  M[idx], mu, phi, Sigma_inv,
                                  max_iter = max_iter, gtol = gtol)
    lambda_hat[, j] <- res$lambda_hat
    S_hat[[j]] <- res$S_hat
    lp_total <- lp_total + res$lp_mode
    ld_S_total <- ld_S_total + res$log_det_shat
    n_iter_vec[j] <- res$n_iter
    grad_norm_vec[j] <- res$grad_norm
  }
  list(lambda_hat = lambda_hat, S_hat = S_hat,
       lp_total = lp_total, ld_S_total = ld_S_total,
       log_det_Sigma = ldSigma,
       n_iter_vec = n_iter_vec, grad_norm_vec = grad_norm_vec)
}

fit_pfa_dense_traced <- function(Y, X, group, K,
                                 M = rowSums(Y),
                                 max_iter = 60, tol = 1e-4,
                                 lambda_phi = 0, sigma2_init = 0.3,
                                 verbose = FALSE,
                                 estep_max_iter = 100, estep_gtol = 1e-3,
                                 B_true = NULL) {
  
  N <- nrow(Y)
  Q <- ncol(Y)
  P <- ncol(X)
  J <- max(group)
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
  es <- NULL
  
  for (em in 1:max_iter) {
    
    es <- estep_dense(J, group, Y, X, M, mu, phi, B, sigma2,
                      max_iter = estep_max_iter, gtol = estep_gtol)
    lambda_hat <- es$lambda_hat
    S_hat <- es$S_hat
    
    log_ev[em] <- es$lp_total - 0.5 * J * es$log_det_Sigma +
      0.5 * es$ld_S_total
    
    mp <- mstep_phi_dense(Y, X, group, lambda_hat, mu, phi,
                          lambda_phi = lambda_phi)
    mu <- mp$mu; phi <- mp$phi
    
    S_obs <- tcrossprod(lambda_hat) / J
    for (j in 1:J) S_obs <- S_obs + S_hat[[j]] / J
    rt <- rubin_thayer_spherical(S_obs, K,
                                 B_init = B, sigma2_init = sigma2)
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
        abs(log_ev[em] - log_ev[em - 1]) < tol) { converged <- TRUE; break }
  }
  
  list(mu = mu, phi = phi, B = B, sigma2 = sigma2,
       lambda_hat = lambda_hat, S_hat = S_hat,
       log_evidence = log_ev[1:em], converged = converged, iterations = em,
       sigma2_trace = sigma2_trace[1:em],
       B_dist_trace = if (!is.null(B_true)) B_dist_trace[1:em] else NULL,
       inner_iter_mean_final = mean(es$n_iter_vec, na.rm = TRUE),
       grad_norm_mean_final = mean(es$grad_norm_vec, na.rm = TRUE),
       inner_maxiter_hit_pct_final = mean(es$n_iter_vec >= estep_max_iter, na.rm = TRUE))
}

all_perms <- function(K) {
  if (K == 1) return(matrix(1, 1, 1))
  sub <- all_perms(K - 1)
  do.call(rbind, lapply(1:K, function(i) cbind(i, sub + (sub >= i))))
}

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

converged_ok <- function(fit, ll_tol, sigma2_tol) {
  le <- fit$log_evidence
  s2 <- fit$sigma2_trace
  n <- length(le)
  if (n < 5) return(FALSE)
  ll_ok <- abs(le[n] - le[n - 1]) < ll_tol
  if (is.null(s2)) return(ll_ok)
  m <- length(s2)
  s2_ok <- abs(s2[m] - s2[m - 1]) < sigma2_tol
  ll_ok && s2_ok
}

first_convergence_iter <- function(fit, ll_tol, sigma2_tol, min_iter = 5) {
  le <- fit$log_evidence
  s2 <- fit$sigma2_trace
  n <- length(le)
  if (n < min_iter)
    return(NA)
  for (i in min_iter:n) {
    ll_ok <- abs(le[i] - le[i - 1]) < ll_tol
    s2_ok <- if (is.null(s2)) TRUE else abs(s2[i] - s2[i - 1]) < sigma2_tol
    if (ll_ok && s2_ok)
      return(i)
  }
  NA
}

run_one_method <- function(method, dat, max_iter, tol) {
  t0 <- proc.time()[3]
  fit <- switch(method,
                base = fit_pfa_dense_traced(dat$Y, dat$X, dat$group, K = dat$K, M = dat$M,
                                            max_iter = max_iter, tol = tol, B_true = dat$true$B),
                normal_step = fit_pfa_woodbury_traced(dat$Y, dat$X, dat$group, K = dat$K, M = dat$M,
                                                      max_iter = max_iter, tol = tol,
                                                      exact_Shat = FALSE, armijo_bt = FALSE, B_true = dat$true$B),
                armijo = fit_pfa_woodbury_traced(dat$Y, dat$X, dat$group, K = dat$K, M = dat$M,
                                                 max_iter = max_iter, tol = tol,
                                                 exact_Shat = FALSE, armijo_bt = TRUE, B_true = dat$true$B)
  )
  runtime <- proc.time()[3] - t0
  tc <- tucker_congruence(fit$B, dat$true$B)
  rv <- rv_coefficient(fit$B, dat$true$B)
  summary_row <- data.frame(
    method = method,
    converged = fit$converged,
    iterations = fit$iterations,
    runtime_sec = as.numeric(runtime),
    log_ev_final = tail(fit$log_evidence, 1),
    sigma2_final = fit$sigma2,
    sigma2_true = dat$true$sigma2,
    sigma2_rel_err = abs(fit$sigma2 - dat$true$sigma2) / dat$true$sigma2,
    B_dist_final = tail(fit$B_dist_trace, 1),
    B_dist_min = min(fit$B_dist_trace),
    tucker_mean = tc$phi_mean,
    rv_coefficient = rv,
    inner_iter_mean_final = fit$inner_iter_mean_final,
    grad_norm_mean_final = fit$grad_norm_mean_final,
    inner_maxiter_hit_pct_final = fit$inner_maxiter_hit_pct_final
  )
  list(summary = summary_row, fit = fit)
}

Q <- c(50)
seeds <- 3
J <- c(50)
K <- 2
N_per_grp <- 15
P <- 3
sigma2 <- 0.3
M_rate <- 150
max_iter <- 100
ll_tol <- 1e-2
sigma2_tol <- 1e-6

methods <- c("base", "normal_step", "armijo")
all_results <- list()

for (q in Q) {
  for (j in J) {
    for (seed in seeds) {
      
      dat <- simulate_pfa_data(Q = q, K = K, J = j, N_per_group = N_per_grp,
                               P = P, sigma2 = sigma2, M_rate = M_rate, seed = seed)
      
      for (method in methods) {
        out <- run_one_method(method, dat, max_iter, tol = 1e-10)
        res <- out$summary
        res$converged_at_max_iter <- converged_ok(out$fit, ll_tol, sigma2_tol)
        res$first_convergence_iter <- first_convergence_iter(out$fit, ll_tol, sigma2_tol)
        res$Q <- q
        res$J <- j
        res$seed <- seed
        
        cat(sprintf("Q=%d J=%d seed=%d method=%s\n", q, j, seed, method))
        print(res)
        
        all_results[[length(all_results) + 1]] <- res
      }
    }
  }
}

results_all <- do.call(rbind, all_results)
rownames(results_all) <- NULL
print(results_all)

# Monotonically check: base
le <- fit_base$log_evidence
d_le <- diff(le)
s2 <- fit_base$sigma2_trace
bd <- fit_base$B_dist_trace

cat("log_ev range:", range(le), "\n")
cat("number of decreasing steps:", sum(d_le < 0), "out of", length(d_le), "\n")
cat("monotone increasing overall?", all(d_le >= 0), "\n")
cat("iter of log_ev peak:", which.max(le), " value:", max(le), "\n")
cat("iter of B_dist min:", which.min(bd), " value:", min(bd), "\n")
cat("log_ev at iter 100:", le[100], " vs peak:", max(le), "\n")
cat("sigma2 at iter 100:", s2[100], " true:", sigma2, "\n")
cat("B_dist at iter 100:", bd[100], " vs min:", min(bd), "\n")

print(round(le, 2))
print(round(s2, 4))
print(round(bd, 4))