library(glmmTMB)
library(bench)

stopifnot(packageVersion("glmmTMB") >= "1.1.8")

source("pfa_base_multinomial.R") # base multinomial
source("pfa_woodbury_only.R") # woodbury
source("pfa_woodbury_lam_corrected.R") # woodbury + lambda correction
source("sim_data.R")
source("basic_functions.R")
# source("test_05_corr_monotone.R") # decouple + corr
# Already pasted in this file
source("pfa_decorr_squarem.R")

fit_glmmTMB_ref <- function(Y, X, group, K, time_limit = Inf, verbose = FALSE) {
  long <- build_long(Y, X, group)
  t0 <- proc.time()["elapsed"]
  
  n_obs_levels <- length(unique(long$obs))
  n_gc_pairs <- length(unique(interaction(long$group, long$category, drop = TRUE)))
  if (n_obs_levels != n_gc_pairs) {
    stop(sprintf("obs is not unique per (group, category): %d obs levels vs %d (group, category) pairs.",
                 n_obs_levels, n_gc_pairs))
  }
  if (anyNA(long$count) || anyNA(long$category) || anyNA(long$group) ||
      anyNA(long$obs) || anyNA(long$log_total)) stop("NA found in build_long() output.")
  if (any(!is.finite(long$log_total))) stop("long$log_total has non-finite value(s).")
  
  form <- as.formula(paste0("count ~ category + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))
  fit_once <- function(ctrl) {
    tryCatch(glmmTMB(form, data = long, family = poisson(link = "log"),
                     offset = long$log_total, control = ctrl),
             error = function(e) e)
  }
  ctrl_default <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  fit <- fit_once(ctrl_default)
  if (inherits(fit, "error")) {
    ctrl_res <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000),
                               start_method = list(method = "res"))
    fit2 <- fit_once(ctrl_res)
    if (inherits(fit2, "error")) {
      el <- proc.time()["elapsed"] - t0
      return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE, converged = FALSE,
                  error = paste0("default: ", conditionMessage(fit), " | res-start: ", conditionMessage(fit2))))
    }
    fit <- fit2
  }
  el <- proc.time()["elapsed"] - t0
  conv <- isTRUE(fit$sdr$pdHess)
  
  L <- tryCatch({
    vc <- as.matrix(VarCorr(fit)$cond$group)
    e  <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
    e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(e$values[1:K], 0)), K)
  }, error = function(e) NULL)
  if (is.null(L)) {
    L <- tryCatch(as.matrix(fit$obj$env$report(fit$fit$parfull)$fact_load[[1]]),
                  error = function(e) NULL)
  }
  if (is.null(L)) {
    return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE, converged = conv,
                error = "loading extraction failed"))
  }
  sigma2_hat <- tryCatch(as.numeric(VarCorr(fit)$cond$obs)[1], error = function(e) NA_real_)
  list(B = L, sigma2 = sigma2_hat, time = el, ok = TRUE, converged = conv, error = NA_character_)
}

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


# fit_pfa_woodbury_G_wrapper <- function(Y, X, group, K, M, max_iter = 80, fix_sigma2 = NULL) {
#   fit_pfa_woodbury_G(Y, X, group, K, M = M, max_iter = max_iter, tol = 1e-8,
#                       sigma2_init = 0.3, verbose = FALSE, fix_sigma2 = fix_sigma2)
# }

# Monotone Check
is_monotone_loglik <- function(fit) {
  all(diff(fit$log_evidence) >= 0)
}

monotone_diagnostics <- function(fit) {
  d <- diff(fit$log_evidence)
  list(monotone = all(d >= 0), n_decreasing = sum(d < 0),
       first_decrease_iter = if (any(d < 0)) which(d < 0)[1] + 1 else NA_integer_)
}

# RUN CODE
Js    <- c(50, 100)
Qs    <- c(50, 100)
# Qs <- c(50)
seeds <- c(1,3) # Observed severe unstable for old woodbury + PME
K <- 2
Nj <- 15

q_results <- data.frame(
  J = integer(), Q = integer(), seed = integer(), method = character(),
  time_s = numeric(), sigma2 = numeric(), d_true = numeric(), converged = logical(),
  log_lik_monotone = logical(),
  stringsAsFactors = FALSE
)

cat(sprintf("%5s %6s %6s | %14s %10s %10s %10s %8s\n %10s\n",
            "J", "Q", "seed", "method", "time(s)", "sigma2", "d vs true", "conv?", "LL mono?"))

for (J in Js) {
  for (Q in Qs) {
    for (seed in seeds) {
      
      # target = (J - K) / J # Large sigma2 test
      target = 0.3 # Small sigma2 test

      dat <- simulate_pfa_data(Q = Q, K = K, J = J, P = 3, N_per_group = Nj, sigma2 = target, M_rate = 150, seed = seed)
      
      # b_base <- bench::mark(
      #   fw_base <<- fit_pfa_dense(dat$Y, dat$X, dat$group, K = K, M = dat$M,
      #                             max_iter = 80, tol = 1e-3, sigma2_init = target, 
      #                             verbose = FALSE, estep_max_iter = 100,
      #                             estep_gtol = 1e-3),
      #   iterations = 1, check = FALSE, memory = FALSE)
      # t_base <-as.numeric(b_base$median)

      b_wb <- bench::mark(
        fw_only <<- fit_pfa_wbonly_traced(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                                     max_iter = 80, tol = 1e-3, sigma2_init = target, 
                                     verbose = FALSE, estep_max_iter = 100,
                                     estep_gtol = 1e-3, B_true = NULL),
        iterations = 1, check = FALSE, memory = FALSE)
      t_wb <- as.numeric(b_wb$median)

      # REML-corrected sigma2
      # s2_reml <- fw_old$sigma2 * J / (J - K)

      # corrected (Woodbury + Edgeworth)
      b_corr <- bench::mark(
        fw_corr <<- fit_pfa_woodbury_lam_corr(dat$Y, dat$X, dat$group, K, M = dat$M,
                                               max_iter = 80, tol = 1e-3, sigma2_init = target,
                                               verbose = FALSE, estep_max_iter = 100, estep_gtol = 1e-3,
                                               B_true = NULL,
                                               corr_max_rel = 0.5, 
                                               use_lambda_correction = TRUE),
        iterations = 1, check = FALSE, memory = FALSE)
      t_corr <- as.numeric(b_corr$median)
      
      b_decorr <- bench::mark(
        fw_decorr <<- fit_pfa_woodbury_lam_corr_decoupled(dat$Y, dat$X, dat$group, K,
                                                          M = dat$M,
                                                          max_iter = 80, tol = 1e-3,
                                                          lambda_phi = 0, sigma2_init = target,
                                                          verbose = FALSE,
                                                          estep_max_iter = 100, estep_gtol = 1e-3,
                                                          B_true = NULL,
                                                          use_lambda_correction = TRUE,
                                                          decouple_muphi_lambda = TRUE,
                                                          corr_max_rel = 0.5,
                                                          fix_sigma2 = NULL,
                                                          freeze_tol = 1e-5, freeze_patience = 3,
                                                          use_aitken = FALSE, aitken_window = 3, aitken_tol = 1e-4),
        iterations = 1, check = FALSE, memory = FALSE)
      t_decorr <- as.numeric(b_decorr$median)
      
      b_decorr_squarem <- bench::mark(
        fw_decorr_squarem <<- fit_pfa_decorr_squarem(dat$Y, dat$X, dat$group, K,
                                                          M = dat$M,
                                                          max_iter = 80, tol = 1e-3,
                                                          lambda_phi = 0, sigma2_init = target,
                                                          verbose = FALSE,
                                                          estep_max_iter = 100, estep_gtol = 1e-3,
                                                          B_true = NULL,
                                                          use_lambda_correction = TRUE,
                                                          decouple_muphi_lambda = TRUE,
                                                          corr_max_rel = 0.5,
                                                          use_squarem = TRUE,
                                                          objfn_inc = 0),
        iterations = 1, check = FALSE, memory = FALSE)
      t_decorr_squarem <- as.numeric(b_decorr_squarem$median)

      # # new (G, unfrozen)
      # b_G <- bench::mark(
      #   fw_G <<- fit_pfa_woodbury_G_wrapper(dat$Y, dat$X, dat$group, K, dat$M,
      #                                        max_iter = 80, fix_sigma2 = NULL),
      #   iterations = 1, check = FALSE, memory = FALSE)
      # t_G <- as.numeric(b_G$median)
      # 
      # # new (G, frozen at REML value)
      # b_G_frz <- bench::mark(
      #   fw_G_frz <<- fit_pfa_woodbury_G_wrapper(dat$Y, dat$X, dat$group, K, dat$M,
      #                                            max_iter = 80, fix_sigma2 = s2_reml),
      #   iterations = 1, check = FALSE, memory = FALSE)
      # t_G_frz <- as.numeric(b_G_frz$median)

      # # glmmTMB rr()
      # b_gt <- tryCatch(
      #   bench::mark(fr <<- fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K = K, verbose = FALSE),
      #               iterations = 1, check = FALSE, memory = FALSE),
      #   error = function(e) NULL)
      # ok <- !is.null(b_gt) && isTRUE(fr$ok)
      # t_gt <- if (ok) as.numeric(b_gt$median) else NA
      
      # d_base_true <- subspace_dist(fw_base$B, dat$true$B)
      d_only_true  <- subspace_dist(fw_only$B, dat$true$B)
      d_corr_true <- subspace_dist(fw_corr$B, dat$true$B)
      d_decorr_true <- subspace_dist(fw_decorr$B, dat$true$B)
      d_decorr_squarem_true <- subspace_dist(fw_decorr_squarem$B, dat$true$B)
      # d_G_true    <- subspace_dist(fw_G$B, dat$true$B)
      # d_Gf_true   <- subspace_dist(fw_G_frz$B, dat$true$B)
      # d_gt_true   <- if (ok) subspace_dist(fr$B, dat$true$B) else NA

      # rows <- data.frame(
      #   J = J, Q = Q, seed = seed,
      #   method = c("old (Woodbury EM)", "corrected (Woodbury+Edgeworth)",
      #              "new (G, unfrozen)", "new (G, frozen)", "glmmTMB rr()"),
      #   time_s = c(t_wb, t_corr, t_G, t_G_frz, t_gt),
      #   d_true = c(d_old_true, d_corr_true, d_G_true, d_Gf_true, d_gt_true),
      #   d_vs_wb_old = c(0,
      #                    subspace_dist(fw_corr$B, fw_old$B),
      #                    subspace_dist(fw_G$B, fw_old$B),
      #                    subspace_dist(fw_G_frz$B, fw_old$B),
      #                    if (ok) subspace_dist(fr$B, fw_old$B) else NA),
      #   converged = c("-", "-", "-", "-", if (ok) ifelse(isTRUE(fr$converged), "yes", "NO") else "FAILED"),
      #   stringsAsFactors = FALSE
      # )
      rows <- data.frame(
        J = J, Q = Q, seed = seed,
        method = c("woodbury only", "corrected", "decorr", "decorr+squarem"),
        time_s = c(t_wb, t_corr, t_decorr, t_decorr_squarem),
        sigma2 = c(fw_only$sigma2, fw_corr$sigma2, fw_decorr$sigma2, fw_decorr_squarem$sigma2),
        d_true = c(d_only_true, d_corr_true, d_decorr_true, d_decorr_squarem_true),
        converged = c(fw_only$converged, fw_corr$converged,
                      fw_decorr$converged, fw_decorr_squarem$converged),
        log_lik_monotone = c(sapply(list(fw_only, fw_corr, fw_decorr, fw_decorr_squarem), function(f) monotone_diagnostics(f)$monotone)),
        stringsAsFactors = FALSE
      )
      q_results <- rbind(q_results, rows)
      
      for (i in seq_len(nrow(rows))) {
        cat(sprintf("%5d %6d %6d | %14s %10.2f %10.4f %10s %8s\n %10s\n",
                    J, Q, seed, rows$method[i], rows$time_s[i], rows$sigma2[i],
                    if (is.na(rows$d_true[i])) "-" else sprintf("%.4f", rows$d_true[i]),
                    if (is.na(rows$converged[i])) "-" else ifelse(rows$converged[i], "yes", "NO"),
                    if (is.na(rows$log_lik_monotone[i])) "-" else ifelse(rows$log_lik_monotone[i], "yes", "NO")))
      }
      cat(sprintf("  J=%d Q=%d seed=%d | true sigma2=%.4f\n", J, Q, seed, dat$true$sigma2))
    }
  }
}

# Calculate the mean of estimations
avg_time  <- aggregate(time_s ~ method + Q + J, data = q_results,
                       FUN = function(x) mean(x, na.rm = TRUE))
avg_sigma <- aggregate(sigma2 ~ method + Q + J, data = q_results,
                       FUN = function(x) mean(x, na.rm = TRUE))
avg_d     <- aggregate(d_true ~ method + Q + J, data = q_results,
                       FUN = function(x) mean(x, na.rm = TRUE))

avg_results <- Reduce(function(a, b) merge(a, b, by = c("method", "Q", "J")),
                      list(avg_time, avg_sigma, avg_d))
avg_results <- avg_results[order(avg_results$Q, avg_results$J, avg_results$method), ]

cat("\n Averaged across seeds, per (method, Q, J) \n")
cat(sprintf("%14s %6s %6s | %10s %10s %10s\n",
            "method", "Q", "J", "time(s)", "sigma2", "d_true"))
for (i in seq_len(nrow(avg_results))) {
  cat(sprintf("%14s %6d %6d | %10.2f %10.4f %10.4f\n",
              avg_results$method[i], avg_results$Q[i], avg_results$J[i],
              avg_results$time_s[i], avg_results$sigma2[i], avg_results$d_true[i]))
}

