cv_select_sigma2 <- function(Y, X, group, K,
                             M = rowSums(Y),
                             sigma2_grid = NULL,   
                             n_folds = 5,
                             fold_seed = 1,
                             cv_max_iter = 40,      
                             cv_tol = 1e-6,
                             refit_max_iter = 150, 
                             refit_tol = 1e-8,
                             probe_max_iter = 150,  
                             probe_tol = 1e-10,      
                             lambda_phi = 0,
                             estep_max_iter = 100, estep_gtol = 1e-3, exact_Shat = FALSE,
                             B_true = NULL,           
                             verbose = TRUE) {
  
  stopifnot(all(X[, 1] == 1))
  J <- max(group)
  
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  
  if (is.null(sigma2_grid)) {
    if (verbose) message("no sigma2_grid supplied -- running auto-detection on the full data first to center a default grid")
    auto_probe <- fit_pfa_woodbury_fixed_sigma2(Y, X, group, K, M = M,
                                                max_iter = probe_max_iter, tol = probe_tol,
                                                fix_sigma2 = "auto", use_aitken = TRUE,
                                                lambda_phi = lambda_phi,
                                                estep_max_iter = estep_max_iter, estep_gtol = estep_gtol,
                                                exact_Shat = exact_Shat, verbose = FALSE)
    if (is.na(auto_probe$sigma2_freeze_value) && verbose) {
      message("WARNING: Aitken never triggered even at probe_max_iter -- falling back to the raw last sigma2 value, which is a weaker estimate. Consider raising probe_max_iter.")
    }
    center <- if (!is.na(auto_probe$sigma2_freeze_value)) auto_probe$sigma2_freeze_value else tail(auto_probe$sigma2_trace, 1)
    sigma2_grid <- center * c(0.6, 0.75, 0.9, 1.0, 1.1, 1.25, 1.4)
    sigma2_grid <- sigma2_grid[sigma2_grid > 0]
    if (verbose) message(sprintf("default grid centered at auto-detected %.4f: %s",
                                 center, paste(sprintf("%.4f", sigma2_grid), collapse = ", ")))
  }
  
  set.seed(fold_seed)
  fold_of_group <- sample(rep(1:n_folds, length.out = J))
  
  results <- data.frame(sigma2 = numeric(), fold = integer(),
                        heldout_ll_per_group = numeric(), ok = logical())
  
  for (cand in sigma2_grid) {
    for (fold in 1:n_folds) {
      
      test_groups  <- sort(which(fold_of_group == fold))
      train_groups <- sort(which(fold_of_group != fold))
      
      train_rows <- group %in% train_groups
      test_rows  <- group %in% test_groups
      
      Y_tr <- Y[train_rows, , drop = FALSE]; X_tr <- X[train_rows, , drop = FALSE]
      M_tr <- M[train_rows]
      group_tr <- as.integer(factor(group[train_rows], levels = train_groups))
      
      Y_te <- Y[test_rows, , drop = FALSE]; X_te <- X[test_rows, , drop = FALSE]
      M_te <- M[test_rows]
      group_te <- as.integer(factor(group[test_rows], levels = test_groups))
      J_te <- length(test_groups)
      
      fit_ok <- TRUE
      heldout_ll <- NA_real_
      
      fit_tr <- tryCatch(
        fit_pfa_woodbury_fixed_sigma2(Y_tr, X_tr, group_tr, K, M = M_tr,
                                      max_iter = cv_max_iter, tol = cv_tol,
                                      lambda_phi = lambda_phi, fix_sigma2 = cand,
                                      estep_max_iter = estep_max_iter, estep_gtol = estep_gtol,
                                      exact_Shat = exact_Shat, verbose = FALSE),
        error = function(e) { fit_ok <<- FALSE; NULL })
      
      if (fit_ok && !is.null(fit_tr)) {
        es_te <- tryCatch(
          estep_wb(J_te, group_te, Y_te, X_te, M_te, fit_tr$mu, fit_tr$phi, fit_tr$B, cand,
                   max_iter = estep_max_iter, gtol = estep_gtol, exact_Shat = exact_Shat),
          error = function(e) { fit_ok <<- FALSE; NULL })
        
        if (fit_ok && !is.null(es_te)) {
          heldout_ll_total <- es_te$lp_total - 0.5 * J_te * es_te$log_det_Sigma + 0.5 * es_te$ld_S_total
          heldout_ll <- heldout_ll_total / J_te   # per held-out group, so folds of unequal size are comparable
        }
      }
      
      results <- rbind(results, data.frame(sigma2 = cand, fold = fold,
                                           heldout_ll_per_group = heldout_ll, ok = fit_ok))
      
      if (verbose) {
        cat(sprintf("  sigma2=%.4f fold=%d/%d: heldout_ll_per_group=%s\n",
                    cand, fold, n_folds, if (fit_ok) sprintf("%.4f", heldout_ll) else "FAILED"))
      }
    }
  }
  
  ok_results <- results[results$ok & !is.na(results$heldout_ll_per_group), ]
  if (nrow(ok_results) == 0) stop("all candidate/fold combinations failed -- check sigma2_grid and data")
  
  mean_ll <- aggregate(heldout_ll_per_group ~ sigma2, data = ok_results, FUN = mean)
  n_ok    <- aggregate(ok ~ sigma2, data = results, FUN = sum)
  summary_tbl <- merge(mean_ll, n_ok, by = "sigma2")
  names(summary_tbl)[names(summary_tbl) == "ok"] <- "n_folds_ok"
  summary_tbl <- summary_tbl[order(summary_tbl$sigma2), ]
  
  best_sigma2 <- summary_tbl$sigma2[which.max(summary_tbl$heldout_ll_per_group)]
  
  if (verbose) {
    cat("\n--- CV summary (mean held-out log-lik per group, higher is better) ---\n")
    print(summary_tbl)
    cat(sprintf("\nselected sigma2 = %.4f\n", best_sigma2))
    cat("refitting on the full data at the selected sigma2...\n")
  }
  
  final_fit <- fit_pfa_woodbury_fixed_sigma2(Y, X, group, K, M = M,
                                             max_iter = refit_max_iter, tol = refit_tol,
                                             lambda_phi = lambda_phi, fix_sigma2 = best_sigma2,
                                             estep_max_iter = estep_max_iter, estep_gtol = estep_gtol,
                                             exact_Shat = exact_Shat, B_true = B_true,
                                             verbose = FALSE)
  
  list(sigma2_selected = best_sigma2,
       cv_table = results,
       cv_summary = summary_tbl,
       fold_of_group = fold_of_group,
       final_fit = final_fit)
}