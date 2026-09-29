source("setup.R")
source("basic_functions.R")

requireNamespace("glmmTMB", quietly = TRUE)
requireNamespace("gllvm", quietly = TRUE)
requireNamespace("COAP", quietly = TRUE)

fit_glmmTMB_ref <- function(Y, X, group, K, n_restart = 5, jitter_sd = 0.2) {
  overall_start <- Sys.time()
  long <- build_long(Y, X, group)

  n_obs_levels <- length(unique(long$obs))
  n_gc_pairs <- length(unique(interaction(long$group, long$category, drop = TRUE)))
  if (n_obs_levels != n_gc_pairs) {
    stop(sprintf("obs is not unique per (group, category): %d obs levels vs %d (group, category) pairs.",
                 n_obs_levels, n_gc_pairs))
  }
  if (anyNA(long$count) || anyNA(long$category) || anyNA(long$group) ||
      anyNA(long$obs) || anyNA(long$log_total)) stop("NA found in build_long() output.")
  if (any(!is.finite(long$log_total))) stop("long$log_total has non-finite value(s).")

  xnames <- paste0("x", 2:ncol(X))
  covar_terms <- paste(paste0("category:", xnames), collapse = " + ")
  form <- as.formula(paste0("count ~ category + ", covar_terms,
                            " + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))
  ctrl <- glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))

  extract <- function(fit, elapsed_sec) {
    conv <- isTRUE(fit$sdr$pdHess)
    L <- tryCatch({
      vc <- as.matrix(glmmTMB::VarCorr(fit)$cond$group)
      e <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
      e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(e$values[1:K], 0)), K)
    }, error = function(e) NULL)
    if (is.null(L)) {
      L <- tryCatch(as.matrix(fit$obj$env$report(fit$fit$parfull)$fact_load[[1]]),
                    error = function(e) NULL)
    }
    if (is.null(L)) {
      return(list(B = NULL, sigma2 = NA_real_, ok = FALSE, converged = conv,
                 loglik = NA_real_, aic = NA_real_, error = "loading extraction failed",
                 elapsed_sec = elapsed_sec))
    }
    L <- apply_PLT(L)
    sigma2_hat <- tryCatch(as.numeric(glmmTMB::VarCorr(fit)$cond$obs)[1], error = function(e) NA_real_)
    loglik_hat <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
    aic_hat <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
    list(B = L, sigma2 = sigma2_hat, ok = TRUE, converged = conv,
        loglik = loglik_hat, aic = aic_hat, error = NA_character_,
        elapsed_sec = elapsed_sec)
  }

  bm0 <- bench::mark(
    fit0 <- tryCatch(glmmTMB::glmmTMB(form, data = long, family = poisson(link = "log"),
                                     offset = log_total, control = ctrl),
                     error = function(e) e),
    memory = FALSE, iterations = 1, check = FALSE)
  fit0_elapsed <- as.numeric(bm0$total_time)

  if (inherits(fit0, "error")) {
    fail <- list(B = NULL, sigma2 = NA_real_, ok = FALSE, converged = FALSE,
                loglik = NA_real_, aic = NA_real_, error = conditionMessage(fit0),
                elapsed_sec = fit0_elapsed)
    all_restarts <- rep(list(fail), n_restart)
    overall_sec <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
    return(list(best = fail, best_idx = NA_integer_, all_restarts = all_restarts,
               timing = list(best_restart_sec = fail$elapsed_sec,
                             all_restarts_sum_sec = sum(sapply(all_restarts, function(r) r$elapsed_sec))),
               meta = list(pbs_job_id = Sys.getenv("PBS_JOBID"), elapsed_sec = overall_sec)))
  }

  restarts <- vector("list", n_restart)
  restarts[[1]] <- extract(fit0, fit0_elapsed)

  if (n_restart > 1) {
    pl0 <- fit0$obj$env$parList()
    rm(fit0)
    gc(verbose = FALSE)
    for (r in 2:n_restart) {
      pl_jit <- pl0
      pl_jit$b <- pl0$b + rnorm(length(pl0$b), 0, jitter_sd)
      pl_jit$beta <- pl0$beta + rnorm(length(pl0$beta), 0, jitter_sd)
      bm_r <- bench::mark(
        fit_r <- tryCatch(glmmTMB::glmmTMB(form, data = long, family = poisson(link = "log"),
                                          offset = log_total, control = ctrl, start = pl_jit),
                          error = function(e) e),
        memory = FALSE, iterations = 1, check = FALSE)
      r_elapsed <- as.numeric(bm_r$total_time)
      restarts[[r]] <- if (inherits(fit_r, "error")) {
        list(B = NULL, sigma2 = NA_real_, ok = FALSE, converged = FALSE,
            loglik = NA_real_, aic = NA_real_, error = conditionMessage(fit_r),
            elapsed_sec = r_elapsed)
      } else extract(fit_r, r_elapsed)
      rm(fit_r)
      gc(verbose = FALSE)
    }
  }

  aics <- sapply(restarts, function(r) if (isTRUE(r$ok) && is.finite(r$aic)) r$aic else Inf)
  best_idx <- if (all(!is.finite(aics))) NA_integer_ else which.min(aics)
  best <- if (is.na(best_idx)) restarts[[1]] else restarts[[best_idx]]

  overall_sec <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
  sum_sec <- sum(sapply(restarts, function(r) r$elapsed_sec))

  list(best = best, best_idx = best_idx, all_restarts = restarts,
      timing = list(best_restart_sec = best$elapsed_sec, all_restarts_sum_sec = sum_sec),
      meta = list(pbs_job_id = Sys.getenv("PBS_JOBID"), elapsed_sec = overall_sec))
}

pick_field <- function(x, candidates) {
  for (nm in candidates) if (!is.null(x[[nm]])) return(x[[nm]])
  NULL
}

orient_B <- function(B, Q, K) {
  B <- as.matrix(B)
  if (nrow(B) != Q && ncol(B) == Q) B <- t(B)
  if (ncol(B) > K) B <- B[, seq_len(K), drop = FALSE]
  B
}

aggregate_to_group <- function(Y, X, group) {
  g <- factor(group)
  nj <- as.vector(table(g))
  list(Y = rowsum(Y, g), X = rowsum(X, g) / nj, n_j = nj)
}

empty_fit <- function(msg, elapsed_sec) {
  list(B = NULL, B_raw = NULL, sigma2 = NA_real_, ok = FALSE, converged = FALSE,
       loglik = NA_real_, aic = NA_real_, u_frac_raw = NA_real_, u_frac_plt = NA_real_,
       error = msg, elapsed_sec = elapsed_sec)
}

pack_fit <- function(L, sigma2, loglik, aic, elapsed_sec, converged, extra = list()) {
  Bp <- apply_PLT(L)
  c(list(B = Bp, B_raw = L, sigma2 = sigma2, ok = TRUE, converged = converged,
         loglik = loglik, aic = aic,
         u_frac_raw = u_frac(L), u_frac_plt = u_frac(Bp),
         error = NA_character_, elapsed_sec = elapsed_sec), extra)
}

wrap_result <- function(restarts, overall_start, extra = list()) {
  aics <- sapply(restarts, function(r) if (isTRUE(r$ok) && is.finite(r$aic)) r$aic else Inf)
  best_idx <- if (all(!is.finite(aics))) NA_integer_ else which.min(aics)
  best <- if (is.na(best_idx)) restarts[[1]] else restarts[[best_idx]]
  list(best = best, best_idx = best_idx, all_restarts = restarts,
       timing = list(best_restart_sec = best$elapsed_sec,
                     all_restarts_sum_sec = sum(sapply(restarts, function(r) r$elapsed_sec))),
       meta = c(list(pbs_job_id = Sys.getenv("PBS_JOBID"),
                     elapsed_sec = as.numeric(difftime(Sys.time(), overall_start, units = "secs"))),
                extra))
}

fit_gllvm_ref <- function(Y, X, group, K, method = "LA",
                          row_eff = c("offset", "fixed"), max_iter = 6000) {
  row_eff <- match.arg(row_eff)
  Yc <- Y
  colnames(Yc) <- paste0("c", seq_len(ncol(Yc)))
  N <- nrow(Yc); Q <- ncol(Yc)
  M <- rowSums(Yc)
  Xd <- as.data.frame(X[, -1, drop = FALSE])
  colnames(Xd) <- paste0("x", 1 + seq_len(ncol(Xd)))
  form <- as.formula(paste("~", paste(colnames(Xd), collapse = " + ")))
  
  args <- list(y = Yc, X = Xd, formula = form, family = poisson(link = "log"),
               num.lv = K, method = method, sd.errors = FALSE, scale.X = FALSE,
               studyDesign = data.frame(group = factor(group)),
               lvCor = ~(1 | group),
               control = list(max.iter = max_iter, maxit = max_iter))
  if (row_eff == "fixed") args$row.eff <- "fixed"
  else args$offset <- matrix(log(M), N, Q)
  
  bm <- bench::mark(
    fit <- tryCatch(do.call(gllvm::gllvm, args), error = function(e) e),
    memory = FALSE, iterations = 1, check = FALSE)
  elapsed <- as.numeric(bm$total_time)
  if (inherits(fit, "error")) return(fail_fit(conditionMessage(fit), elapsed))
  
  L <- tryCatch(as.matrix(gllvm::getLoadings(fit)), error = function(e) NULL)
  if (is.null(L)) {
    L <- tryCatch({
      th <- as.matrix(fit$params$theta)[, seq_len(K), drop = FALSE]
      sv <- fit$params$sigma.lv
      if (is.null(sv)) th else th %*% diag(sv[seq_len(K)], K)
    }, error = function(e) NULL)
  }
  if (is.null(L)) return(fail_fit("loading extraction failed", elapsed))
  L <- as.matrix(L)[, seq_len(K), drop = FALSE]
  if (nrow(L) != Q) L <- t(L)
  
  b0 <- unname(fit$params$beta0)
  phi <- rbind(0, t(unname(as.matrix(fit$params$Xcoef))))
  ll <- as.numeric(fit$logL)
  nlv <- tryCatch(nrow(as.matrix(gllvm::getLV(fit))), error = function(e) NA_integer_)
  
  list(B = apply_PLT(L), B_raw = L, sigma2 = NA_real_,
       mu = b0 - mean(b0), phi = phi,
       omega_hat = eigen(crossprod(apply_PLT(L)), symmetric = TRUE, only.values = TRUE)$values / Q,
       u_frac_raw = u_frac(L), colnorm_raw = sqrt(colSums(L^2)),
       n_lv_rows = nlv, sigma_lv = fit$params$sigma.lv,
       ok = TRUE, converged = is.finite(ll), error = NA_character_, loglik = ll,
       meta = list(pbs_job_id = Sys.getenv("PBS_JOBID"), elapsed_sec = elapsed,
                   row_eff = row_eff, method = method))
}
# COAP

fit_coap_ref <- function(Y, X, group, K, level = c("group", "obs"),
                         rank_use = NULL, n_offset_iter = 5,
                         maxIter = 100, epsELBO = 1e-6) {
  overall_start <- Sys.time()
  level <- match.arg(level)
  
  if (level == "group") {
    g <- factor(group)
    nj <- as.vector(table(g))
    Yc <- rowsum(Y, g)
    Zc <- rowsum(X, g) / nj
  } else {
    nj <- NULL
    Yc <- Y
    Zc <- X
  }
  n <- nrow(Yc)
  Q <- ncol(Yc)
  P <- ncol(Zc)
  M <- rowSums(Yc)
  if (is.null(rank_use)) rank_use <- P
  if (n <= Q) warning(sprintf("COAP fitted with n = %d rows and Q = %d columns.", n, Q))
  
  fac <- M
  fit <- NULL
  el_total <- 0
  offset_trace <- numeric(0)
  
  for (it in seq_len(max(1L, n_offset_iter))) {
    bm <- bench::mark(
      f_it <- tryCatch(COAP::RR_COAP(X_count = Yc, multiFac = fac, Z = Zc,
                                     rank_use = rank_use, q = K,
                                     epsELBO = epsELBO, maxIter = maxIter,
                                     verbose = FALSE),
                       error = function(e) e),
      memory = FALSE, iterations = 1, check = FALSE)
    el_total <- el_total + as.numeric(bm$total_time)
    if (inherits(f_it, "error")) {
      return(wrap_result(list(empty_fit(conditionMessage(f_it), el_total)), overall_start,
                         list(level = level, n_offset_done = it - 1L)))
    }
    fit <- f_it
    if (it == n_offset_iter) break
    Bh <- orient_B(fit$B, Q, K)
    Hh <- as.matrix(fit$H)
    bb <- as.matrix(fit$bbeta)
    if (nrow(bb) != Q) bb <- t(bb)
    eta <- Zc %*% t(bb) + Hh %*% t(Bh)
    fac_new <- M / exp(row_logsumexp(eta))
    offset_trace <- c(offset_trace, mean(abs(log(fac_new) - log(fac))))
    fac <- fac_new
  }
  
  L_raw <- orient_B(fit$B, Q, K)
  L <- apply_PLT(L_raw)
  Hh <- as.matrix(fit$H)
  bb <- as.matrix(fit$bbeta)
  if (nrow(bb) != Q) bb <- t(bb)
  
  mu_hat <- bb[, 1] - mean(bb[, 1])
  phi_hat <- rbind(0, t(bb[, -1, drop = FALSE]))
  
  s2 <- NA_real_
  if (!is.null(fit$invLambda)) s2 <- mean(1 / as.numeric(fit$invLambda))
  else if (!is.null(fit$Lambda)) s2 <- mean(as.numeric(fit$Lambda))
  
  nit <- length(fit$ELBO_seq)
  elbo <- if (nit > 0) as.numeric(tail(fit$ELBO_seq, 1)) else as.numeric(fit$ELBO)
  
  omega_hat <- eigen(crossprod(L), symmetric = TRUE, only.values = TRUE)$values / Q
  
  res <- list(B = L, B_raw = L_raw, H = Hh, bbeta = bb,
              mu = mu_hat, phi = phi_hat, sigma2 = s2,
              ok = TRUE, converged = nit < maxIter,
              loglik = elbo, aic = -2 * elbo,
              u_frac_raw = u_frac(L_raw), u_frac_plt = u_frac(L),
              omega_hat = omega_hat,
              colnorm_raw = sqrt(colSums(L_raw^2)),
              HtH_over_n = crossprod(Hh) / nrow(Hh),
              n_iter = nit, hit_cap = nit >= maxIter,
              n_bad_invLambda = sum(!is.finite(as.numeric(fit$invLambda))),
              eta_range = range(Zc %*% t(bb) + Hh %*% t(L_raw)),
              error = NA_character_, elapsed_sec = el_total)
  
  wrap_result(list(res), overall_start,
              list(level = level, rank_use = rank_use, n_rows = n,
                   n_j = nj, offset_trace = offset_trace))
}