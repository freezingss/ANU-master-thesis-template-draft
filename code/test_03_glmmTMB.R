library(glmmTMB) # for rr() comparison on the level of Poisson surrogate
library(bench) # timing
library(kableExtra)

stopifnot(packageVersion("glmmTMB") >= "1.1.8")   # covariance rr() requires

source("pfa_woodbury.R")
source("sim_data.R")
source("basic_functions.R")

# fit glmmTMB reduced-rank reference
fit_glmmTMB_ref <- function(Y, X, group, K, time_limit = Inf, verbose = TRUE) {
  
  long <- build_long(Y, X, group)
  t0 <- proc.time()["elapsed"]
  
  n_obs_levels <- length(unique(long$obs))
  n_gc_pairs   <- length(unique(interaction(long$group, long$category,
                                            drop = TRUE)))
  if (n_obs_levels != n_gc_pairs) {
    stop(sprintf(
      "obs is not unique per (group, category): %d obs levels vs %d (group, category) pairs. Fix build_long() (e.g. obs <- interaction(group, category, drop = TRUE)).",
      n_obs_levels, n_gc_pairs))
  }

  if (anyNA(long$count) || anyNA(long$category) || anyNA(long$group) ||
      anyNA(long$obs) || anyNA(long$log_total)) {
    stop("NA found in long$count / category / group / obs / log_total - check build_long().")
  }
  if (any(!is.finite(long$log_total))) {
    stop(sprintf(
      "long$log_total has %d non-finite value(s) (likely log(0) from a zero total count). Check for empty groups/categories in Y.",
      sum(!is.finite(long$log_total))))
  }
  
  form <- as.formula(paste0(
    "count ~ category + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))
  
  fit_once <- function(ctrl) {
    tryCatch(
      glmmTMB(form, data = long, family = poisson(link = "log"),
              offset = long$log_total, control = ctrl),
      error = function(e) e)  
  }
  
  ctrl_default <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  fit <- fit_once(ctrl_default)
  
  if (inherits(fit, "error")) {
    msg1 <- conditionMessage(fit)
    if (verbose) message("glmmTMB (default start) failed: ", msg1,
                         " -- retrying with start_method = 'res'")
    
    ctrl_res <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000),
                               start_method = list(method = "res"))
    fit2 <- fit_once(ctrl_res)
    
    if (inherits(fit2, "error")) {
      msg2 <- conditionMessage(fit2)
      if (verbose) message("glmmTMB (start_method='res') also failed: ", msg2)
      el <- proc.time()["elapsed"] - t0
      return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE,
                  converged = FALSE,
                  error = paste0("default: ", msg1, " | res-start: ", msg2)))
    }
    fit <- fit2
  }
  
  el <- proc.time()["elapsed"] - t0
  conv <- isTRUE(fit$sdr$pdHess)
  if (verbose && !conv) {
    message("glmmTMB fit completed but pdHess = FALSE (Hessian not positive definite) -- treat this fit's estimates with caution.")
  }

  L <- tryCatch({
    vc <- VarCorr(fit)$cond$group
    if (is.null(vc)) stop("VarCorr(fit)$cond$group is NULL")
    vc <- as.matrix(vc)
    e  <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
    e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(e$values[1:K], 0)), K)
  }, error = function(e) NULL)
  
  if (is.null(L)) {
    if (verbose) message("VarCorr(fit)$cond$group route failed, falling back to fact_load report")
    L <- tryCatch({
      rep_obj <- fit$obj$env$report(fit$fit$parfull)
      as.matrix(rep_obj$fact_load[[1]])
    }, error = function(e) NULL)
  }
  
  if (is.null(L)) {
    if (verbose) message("Could not extract loadings from either VarCorr() or fact_load; fit itself succeeded")
    return(list(B = NULL, sigma2 = NA_real_, time = el, ok = FALSE,
                converged = conv, error = "loading extraction failed"))
  }

  sigma2_hat <- tryCatch({
    vc_obs <- VarCorr(fit)$cond$obs
    as.numeric(vc_obs)[1]
  }, error = function(e) NA_real_)
  
  list(B = L, sigma2 = sigma2_hat, time = el, ok = TRUE,
       converged = conv, error = NA_character_)
}

long <- build_long_gT(dat$Y, dat$X, dat$group)
length(unique(long$obs)) == length(unique(interaction(long$group, long$category, drop = TRUE))) # TRUE
any(!is.finite(long$log_total)) # FALSE

cat(sprintf("%5s %5s %5s | %14s %14s %16s\n",
            "Q", "K", "J", "wb vs true", "glmmTMB vs true", "wb vs glmmTMB"))

for (Q in c(20, 30, 50)) {
  K <- 2; J <- 30; Nj <- 15
  dat <- simulate_pfa_data(Q = Q, K = K, J = J, N_per_group = Nj, seed = 11)
  
  fw <- fit_pfa_woodbury(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                         max_iter = 30, verbose = FALSE, exact_Shat = FALSE)
  fr <- fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K = K)
  
  d_wb_true <- subspace_dist(fw$B, dat$true$B)
  if (fr$ok) {
    d_gt_true <- subspace_dist(fr$B, dat$true$B)
    d_wb_gt   <- subspace_dist(fw$B, fr$B)
    cat(sprintf("%5d %5d %5d | %14.4f %14.4f %16.4f\n",
                Q, K, J, d_wb_true, d_gt_true, d_wb_gt))
  } else {
    cat(sprintf("%5d %5d %5d | %14.4f %14s %16s\n",
                Q, K, J, d_wb_true, "FAILED", "-"))
  }
}

cat(sprintf("%6s | %14s %14s %10s %10s %10s %8s\n",
            "Q", "Woodbury (s)", "glmmTMB (s)",
            "wb vs gT", "gT vs true",  "wb vs true", "conv?"))

Qs <- c(20, 30, 50)

results <- data.frame(
  Q = numeric(),
  Woodbury_Time = numeric(),
  Woodbury_Memory = numeric(),
  glmmTMB_Time = numeric(),
  glmmTMB_Memory = numeric(),
  Subspace_Distance = numeric(),
  glmmTMB_vs_True = numeric(),
  Woodbury_vs_True = numeric(),
  Converged = character()
)

for (Q in Qs) {
  K <- 2; J <- 30; Nj <- 15
  dat <- simulate_pfa_data(Q = Q, K = K, J = J, N_per_group = Nj, seed = Q)
  
  b_wb <- bench::mark(
    fw <<- fit_pfa_woodbury(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                            max_iter = 30, verbose = FALSE,
                            exact_Shat = FALSE),
    iterations = 3, check = FALSE, memory = TRUE)
  t_wb <- as.numeric(b_wb$median)
  m_wb <- as.numeric(b_wb$mem_alloc) / 1024^2
  
  b_gt <- tryCatch(
    bench::mark(fr <<- fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K = K),
                iterations = 2, check = FALSE, memory = TRUE),
    error = function(e) NULL)
  ok <- !is.null(b_gt) && isTRUE(fr$ok)
  t_gt <- if (ok) as.numeric(b_gt$median) else NA
  m_gt <- if (ok) as.numeric(b_gt$mem_alloc) / 1024^2 else NA
  
  d_wb_gt   <- if (ok) subspace_dist(fw$B, fr$B) else NA
  d_gt_true <- if (ok) subspace_dist(fr$B, dat$true$B) else NA
  d_wb_true <- subspace_dist(fw$B, dat$true$B)
  conv <- if (ok) isTRUE(fr$converged) else FALSE
  
  cat(sprintf("%6d | %14.2f %14s %10s %10s %10s %8s\n",
              Q, t_wb,
              if (ok) sprintf("%.2f", t_gt) else "FAILED/slow",
              if (is.na(d_wb_gt))   "-" else sprintf("%.4f", d_wb_gt),
              if (is.na(d_gt_true)) "-" else sprintf("%.4f", d_gt_true),
              sprintf("%.4f", d_wb_true),
              if (ok) ifelse(conv, "yes", "NO") else "-"))
  
  results <- rbind(
    results,
    data.frame(
      Q = Q,
      Woodbury_Time = t_wb,
      Woodbury_Memory = m_wb,
      glmmTMB_Time = t_gt,
      glmmTMB_Memory = m_gt,
      Subspace_Distance = d_wb_gt,
      glmmTMB_vs_True = d_gt_true,
      Woodbury_vs_True = d_wb_true,
      Converged = ifelse(ok, ifelse(conv, "yes", "no"), "failed")
    )
  )
}

# More seeds
seeds  <- 101:105        
Qs_chk <- c(30, 50)        
K <- 2; J <- 30; Nj <- 15

seed_results <- data.frame(
  Q = integer(), seed = integer(),
  d_wb_true = numeric(), d_gt_true = numeric(), d_wb_gt = numeric(),
  gt_converged = character(), gt_error = character(),
  stringsAsFactors = FALSE
)

cat(sprintf("%5s %6s | %10s %10s %10s %6s\n",
            "Q", "seed", "wb vs true", "gt vs true", "wb vs gt", "conv?"))

for (Qv in Qs_chk) {
  for (s in seeds) {
    dat <- simulate_pfa_data(Q = Qv, K = K, J = J, N_per_group = Nj, seed = s)
    
    fw <- fit_pfa_woodbury(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                           max_iter = 30, verbose = FALSE, exact_Shat = TRUE)
    fr <- fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K = K, verbose = FALSE)
    
    d_wb_true <- subspace_dist(fw$B, dat$true$B)
    d_gt_true <- if (fr$ok) subspace_dist(fr$B, dat$true$B) else NA
    d_wb_gt   <- if (fr$ok) subspace_dist(fw$B, fr$B)        else NA
    
    cat(sprintf("%5d %6d | %10s %10s %10s %6s\n",
                Qv, s,
                sprintf("%.4f", d_wb_true),
                if (fr$ok) sprintf("%.4f", d_gt_true) else "FAILED",
                if (fr$ok) sprintf("%.4f", d_wb_gt)   else "-",
                if (fr$ok) ifelse(isTRUE(fr$converged), "yes", "NO") else "-"))
    
    if (!fr$ok && !is.null(fr$error) && !is.na(fr$error)) {
      message(sprintf("  Q=%d seed=%d glmmTMB failed: %s", Qv, s, fr$error))
    }
    
    seed_results <- rbind(seed_results, data.frame(
      Q = Qv, seed = s,
      d_wb_true = d_wb_true, d_gt_true = d_gt_true, d_wb_gt = d_wb_gt,
      gt_converged = if (fr$ok) ifelse(isTRUE(fr$converged), "yes", "no") else "failed",
      gt_error = if (fr$ok) NA_character_ else fr$error,
      stringsAsFactors = FALSE
    ))
  }
}

for (Qv in Qs_chk) {
  sub  <- seed_results[seed_results$Q == Qv & !is.na(seed_results$d_gt_true), ]
  diff <- sub$d_wb_true - sub$d_gt_true
  cat(sprintf("Q=%d: n=%d valid seeds, mean diff=%.4f, sd=%.4f, glmmTMB closer in %d/%d\n",
              Qv, length(diff), mean(diff), sd(diff), sum(diff > 0), length(diff)))
}