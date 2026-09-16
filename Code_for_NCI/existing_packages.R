source("setup.R")
source("basic_functions.R")

# glmmTMB rr() factor analysis with plug-in offset
fit_glmmTMB_ref <- function(Y, X, group, K, verbose = FALSE) {
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
  
  # form <- as.formula(paste0("count ~ category + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))
  xnames <- paste0("x", 2:ncol(X))
  covar_terms <- paste(paste0("category:", xnames), collapse = " + ")
  form <- as.formula(paste0("count ~ category + ", covar_terms,
                            " + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))
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
      return(list(B = NULL, sigma2 = NA_real_, ok = FALSE, converged = FALSE,
                  error = paste0("default: ", conditionMessage(fit), " | res-start: ", conditionMessage(fit2))))
    }
    fit <- fit2
  }
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
    return(list(B = NULL, sigma2 = NA_real_, ok = FALSE, converged = conv,
                error = "loading extraction failed"))
  }
  sigma2_hat <- tryCatch(as.numeric(VarCorr(fit)$cond$obs)[1], error = function(e) NA_real_)
  list(B = L, sigma2 = sigma2_hat, ok = TRUE, converged = conv, error = NA_character_)
}

# gllvm (may not include covariates)
