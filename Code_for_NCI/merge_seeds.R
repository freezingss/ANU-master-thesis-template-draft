res_dir <- "results"
out_dir <- "results_merged"
dir.create(out_dir, showWarnings = FALSE)

J_grid <- c(100)
Q_grid <- c(50, 100, 150, 200, 300, 400)
seeds <- 1:100
methods <- c("oracle", "wb")

for (J_sel in J_grid) {
  for (Q_sel in Q_grid) {
    merged <- vector("list", length(methods))
    names(merged) <- methods
    for (m in methods) {
      per_method <- vector("list", length(seeds))
      names(per_method) <- paste0("seed", seeds)
      n_found <- 0
      for (i in seq_along(seeds)) {
        path <- file.path(res_dir, sprintf("J%d_Q%d_seed%d_%s.rds", J_sel, Q_sel, seeds[i], m))
        if (file.exists(path)) {
          r <- tryCatch(readRDS(path), error = function(e) {
            cat(sprintf("CORRUPTED (failed to read): %s -- %s\n", basename(path), conditionMessage(e)))
            NULL
          })
          if (!is.null(r)) {
            per_method[[i]] <- r
            n_found <- n_found + 1
          }
        } else {
          cat(sprintf("missing: %s\n", basename(path)))
        }
      }
      merged[[m]] <- per_method
      cat(sprintf("J=%d Q=%d method=%s: %d/%d seeds found\n", J_sel, Q_sel, m, n_found, length(seeds)))
    }
    out_path <- file.path(out_dir, sprintf("J%d_Q%d_merged.rds", J_sel, Q_sel))
    saveRDS(merged, out_path)
    cat(sprintf("wrote %s\n", out_path))
  }
}
