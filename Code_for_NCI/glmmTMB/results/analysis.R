source("basic_functions.R")

gauge_center <- function(B) sweep(as.matrix(B), 2, colMeans(as.matrix(B)))

u_frac_of <- function(B) {
  B <- as.matrix(B)
  fn <- norm(B, "F")
  if (!is.finite(fn) || fn == 0) return(NA_real_)
  sqrt(nrow(B) * sum(colMeans(B)^2)) / fn
}

list_seed_files_glmm <- function(res_dir, J_sel, Q_sel, seeds, method_sel = "glmmTMB") {
  paths <- file.path(res_dir, sprintf("J%d_Q%d_seed%d_%s.rds", J_sel, Q_sel, seeds, method_sel))
  found <- file.exists(paths)
  list(files = paths[found], seeds = seeds[found])
}

extract_one_glmm <- function(path, seed_val) {
  r <- readRDS(path)
  fit <- r$fit
  truth <- r$truth
  ok <- isTRUE(fit$ok)
  metric_names <- c("d_true", "ang_max", "tucker", "rv", "rmse", "sig_err")
  if (ok && !is.null(fit$B)) {
    m <- B_metrics(fit$B, fit$sigma2, truth$B, truth$sigma2)
    mp <- B_metrics(gauge_center(fit$B), fit$sigma2, truth$B, truth$sigma2)
    uf_hat <- u_frac_of(fit$B)
    uf_true <- u_frac_of(truth$B)
  } else {
    m <- setNames(rep(NA_real_, length(metric_names)), metric_names)
    mp <- m
    uf_hat <- NA_real_
    uf_true <- NA_real_
  }
  row <- data.frame(seed = seed_val, ok = ok,
                    converged = isTRUE(fit$converged),
                    error = if (is.null(fit$error) || is.na(fit$error)) NA_character_ else fit$error,
                    sigma2_hat = if (ok) fit$sigma2 else NA_real_,
                    sigma2_true = truth$sigma2)
  for (nm in metric_names) row[[nm]] <- unname(m[nm])
  row$d_true_proj <- unname(mp["d_true"])
  row$ang_max_proj <- unname(mp["ang_max"])
  row$u_frac_hat <- uf_hat
  row$u_frac_true <- uf_true
  row$d_true_drop_pct <- 100 * (row$d_true - row$d_true_proj) / row$d_true
  row$sigma2_bias_pct <- 100 * (row$sigma2_hat - row$sigma2_true) / row$sigma2_true
  row
}

analyse_JQ_glmmTMB <- function(J_sel, Q_sel, res_dir = ".", out_dir = "analysis",
                               expected_seeds = 1:5, verbose = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tag <- sprintf("J%d_Q%d_glmmTMB", J_sel, Q_sel)
  
  sf <- list_seed_files_glmm(res_dir, J_sel, Q_sel, expected_seeds)
  cat(sprintf("%d files found for %s\n", length(sf$files), tag))
  cat("missing seeds:", setdiff(expected_seeds, sf$seeds), "\n")
  
  if (length(sf$files) == 0) stop(sprintf("no files found for %s in %s", tag, normalizePath(res_dir, mustWork = FALSE)))
  
  seed_tbl <- do.call(rbind, Map(extract_one_glmm, sf$files, sf$seeds))
  rownames(seed_tbl) <- NULL
  
  if (verbose) {
    cat("\nper-seed results\n")
    print(seed_tbl)
  }
  
  n_ok <- sum(seed_tbl$ok)
  n_converged <- sum(seed_tbl$converged, na.rm = TRUE)
  cat(sprintf("ok: %d/%d, converged: %d/%d\n", n_ok, nrow(seed_tbl), n_converged, nrow(seed_tbl)))
  if (any(!seed_tbl$ok)) {
    cat("failed seeds and their errors:\n")
    print(seed_tbl[!seed_tbl$ok, c("seed", "error")])
  }
  
  metric_cols <- c("d_true", "d_true_proj", "d_true_drop_pct",
                   "ang_max", "ang_max_proj", "u_frac_hat", "u_frac_true",
                   "tucker", "rv", "rmse", "sig_err", "sigma2_hat", "sigma2_bias_pct")
  use_rows <- seed_tbl[seed_tbl$ok, , drop = FALSE]
  summary_tbl <- do.call(rbind, lapply(metric_cols, function(cn) {
    v <- use_rows[[cn]]
    n <- sum(is.finite(v))
    data.frame(J = J_sel, Q = Q_sel, method = "glmmTMB", quantity = cn,
               n = n,
               mean = if (n > 0) mean(v, na.rm = TRUE) else NA_real_,
               sd = if (n > 1) sd(v, na.rm = TRUE) else NA_real_,
               median = if (n > 0) median(v, na.rm = TRUE) else NA_real_,
               min = if (n > 0) min(v, na.rm = TRUE) else NA_real_,
               max = if (n > 0) max(v, na.rm = TRUE) else NA_real_)
  }))
  summary_tbl$mcse <- summary_tbl$sd / sqrt(summary_tbl$n)
  
  if (verbose) {
    cat("\nseed-averaged key quantities\n")
    print(summary_tbl, row.names = FALSE)
  }
  
  write.csv(seed_tbl, file.path(out_dir, sprintf("per_seed_%s.csv", tag)), row.names = FALSE)
  write.csv(summary_tbl, file.path(out_dir, sprintf("summary_%s.csv", tag)), row.names = FALSE)
  
  list(seed_tbl = seed_tbl, summary_tbl = summary_tbl)
}

collect_grid_glmm <- function(out_dir, J_grid, Q_grid) {
  jq_run <- expand.grid(J = J_grid, Q = Q_grid)
  files <- file.path(out_dir, sprintf("summary_J%d_Q%d_glmmTMB.csv", jq_run$J, jq_run$Q))
  files <- files[file.exists(files)]
  do.call(rbind, lapply(files, read.csv))
}

res_dir <- "."
out_dir <- file.path(res_dir, "analysis")
J_grid <- c(50)
Q_grid <- c(100)
expected_seeds <- 1:10

jq_run <- expand.grid(J = J_grid, Q = Q_grid)
res_all <- list()
for (i in seq_len(nrow(jq_run))) {
  J_i <- jq_run$J[i]
  Q_i <- jq_run$Q[i]
  tag_i <- sprintf("J%d_Q%d_glmmTMB", J_i, Q_i)
  cat(sprintf("\n==================== %s ====================\n", tag_i))
  res_all[[tag_i]] <- tryCatch(
    analyse_JQ_glmmTMB(J_i, Q_i, res_dir = res_dir, out_dir = out_dir, expected_seeds = expected_seeds),
    error = function(e) {
      message(sprintf("%s failed: %s", tag_i, conditionMessage(e)))
      NULL
    }
  )
}

cat("\nfinished:", names(res_all)[!sapply(res_all, is.null)], "\n")
cat("failed:", names(res_all)[sapply(res_all, is.null)], "\n")

summary_grid <- collect_grid_glmm(out_dir, J_grid, Q_grid)
if (nrow(summary_grid) > 0) {
  summary_grid <- summary_grid[order(summary_grid$J, summary_grid$Q, summary_grid$quantity), ]
  write.csv(summary_grid, file.path(out_dir, "summary_grid_glmmTMB.csv"), row.names = FALSE)
  cat("\nseed-averaged key quantities across the grid\n")
  print(summary_grid, row.names = FALSE)
}

Q_test <- c(50,100,150,200,300,400)
J_test<- c(100,200)
for (q in Q_test){
  for (j in J_test){
    print(q/j)
  }
}

