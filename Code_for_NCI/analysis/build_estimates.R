source("basic_functions.R")

res_dirs <- c("G:/J100_Q_seed10_COAP")
out_dir <- c("G:/J100_Q_seed100_COAP_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

J_grid <- c(100)
Q_grid <- c(50, 100, 150, 200, 300, 400)
seeds <- 1:100
methods <- c("coap_group", "coap_obs")
n_cores <- max(1L, as.integer(Sys.getenv("PBS_NCPUS", "1")))
only_expected <- TRUE

array_fields <- c("B", "B_last", "mu", "phi")
metric_names <- c("d_true", "ang_max", "tucker", "rv", "rmse", "sig_err")

ev_of <- function(B) {
  if (is.null(B)) return(NULL)
  eigen(crossprod(B), symmetric = TRUE, only.values = TRUE)$values
}

coap_field <- function(fit, name, alt = name) {
  if (!is.null(fit$best)) fit$best[[name]] else if (!is.null(fit[[name]])) fit[[name]] else fit[[alt]]
}

null2na <- function(x) if (is.null(x) || length(x) != 1) NA else x

est_views <- list(
  wb = function(fit) list(
    ok = !is.null(fit$B) && all(is.finite(fit$B)),
    converged = isTRUE(fit$converged),
    B = fit$B, sigma2 = fit$sigma2, mu = fit$mu, phi = fit$phi,
    B_last = fit$checkpoints$final$B,
    sigma2_last = fit$checkpoints$final$sigma2,
    obj_value = fit$best_log_ev,
    elapsed_sec = fit$meta$elapsed_sec),
  oracle = function(fit) list(
    ok = fit$ok, converged = fit$converged,
    B = fit$B, sigma2 = fit$sigma2,
    elapsed_sec = fit$meta$elapsed_sec),
  glmmTMB = function(fit) list(
    ok = fit$best$ok, converged = fit$best$converged,
    B = fit$best$B, sigma2 = fit$best$sigma2,
    obj_value = fit$best$loglik, aic = fit$best$aic,
    elapsed_sec = fit$meta$elapsed_sec,
    best_restart_sec = fit$timing$best_restart_sec,
    all_restarts_sum_sec = fit$timing$all_restarts_sum_sec,
    n_restart = length(fit$all_restarts), best_idx = fit$best_idx),
  coap_group = function(fit) list(
    ok = isTRUE(coap_field(fit, "ok")),
    converged = isTRUE(coap_field(fit, "converged")),
    B = coap_field(fit, "B"),
    sigma2 = coap_field(fit, "sigma2"),
    mu = coap_field(fit, "mu"),
    obj_value = coap_field(fit, "obj_value"),
    n_iter = coap_field(fit, "n_iter", alt = "iterations"),
    hit_cap = coap_field(fit, "hit_cap"),
    elapsed_sec = fit$meta$elapsed_sec
  ),
  coap_obs = function(fit) list(
    ok = isTRUE(coap_field(fit, "ok")),
    converged = isTRUE(coap_field(fit, "converged")),
    B = coap_field(fit, "B"),
    sigma2 = coap_field(fit, "sigma2"),
    mu = coap_field(fit, "mu"),
    phi = coap_field(fit, "phi"),
    obj_value = coap_field(fit, "obj_value"),
    n_iter = coap_field(fit, "n_iter", alt = "iterations"),
    hit_cap = coap_field(fit, "hit_cap"),
    elapsed_sec = fit$meta$elapsed_sec
  )
)

default_view <- function(fit) list(
  ok = fit$ok, converged = fit$converged,
  B = fit$B, sigma2 = fit$sigma2, mu = fit$mu, phi = fit$phi,
  elapsed_sec = fit$meta$elapsed_sec)

u_frac_of <- function(B) {
  if (is.null(B)) return(NA_real_)
  sum(colSums(B)^2) / (nrow(B) * sum(B^2))
}

bind_fill <- function(lst) {
  cols <- unique(unlist(lapply(lst, names)))
  lst <- lapply(lst, function(d) {
    for (m in setdiff(cols, names(d))) d[[m]] <- NA
    d[cols]
  })
  do.call(rbind, lst)
}

key_of <- function(d) paste(d$J, d$Q, d$seed, d$method, sep = "|")
key3_of <- function(d) paste(d$J, d$Q, d$seed, sep = "|")

list_result_files <- function(dirs) {
  files <- unlist(lapply(dirs, function(d) list.files(d, pattern = "\\.rds$", full.names = TRUE)))
  b <- basename(files)
  parts <- regmatches(b, regexec("^J([0-9]+)_Q([0-9]+)_seed([0-9]+)_(.+)\\.rds$", b))
  hit <- lengths(parts) == 5
  cat(sprintf("%d .rds files found, %d match the naming pattern, %d ignored\n",
              length(files), sum(hit), sum(!hit)))
  if (!any(hit)) stop("no result files matched the naming pattern")
  m <- do.call(rbind, parts[hit])
  data.frame(path = files[hit], J = as.integer(m[, 2]), Q = as.integer(m[, 3]),
             seed = as.integer(m[, 4]), method = m[, 5], stringsAsFactors = FALSE)
}

read_one <- function(path, J, Q, seed, method) {
  out <- list(row = data.frame(J = J, Q = Q, seed = seed, method = method,
                               file_exists = TRUE, stringsAsFactors = FALSE),
              arrays = list(), truth = NULL)
  r <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(r, "error")) {
    out$row$read_ok <- FALSE
    out$row$read_error <- conditionMessage(r)
    return(out)
  }
  out$row$read_ok <- TRUE
  out$row$read_error <- NA_character_
  K <- null2na(r$config$K)
  out$row$K <- K
  view_fn <- if (is.null(est_views[[method]])) default_view else est_views[[method]]
  view <- tryCatch(view_fn(r$fit), error = function(e) e)
  if (inherits(view, "error")) {
    out$row$read_error <- paste("view failed:", conditionMessage(view))
    view <- list()
  }
  for (nm in setdiff(names(view), array_fields)) out$row[[nm]] <- null2na(view[[nm]])
  out$arrays <- view[intersect(names(view), array_fields)]
  B <- view$B
  out$row$B_dim_ok <- isTRUE(!is.null(B) && nrow(B) == Q && ncol(B) == K)
  out$row$u_frac <- u_frac_of(B)
  Bt <- r$truth$B
  s2t <- null2na(r$truth$sigma2)
  s2h <- null2na(view$sigma2)
  m <- setNames(rep(NA_real_, length(metric_names)), metric_names)
  if (!is.null(B) && !is.null(Bt)) {
    m <- tryCatch(B_metrics(B, s2h, Bt, s2t), error = function(e) m)
  }
  for (nm in metric_names) out$row[[nm]] <- unname(m[nm])
  out$row$sigma2_bias_pct <- 100 * (s2h - s2t) / s2t
  ev <- ev_of(B)
  for (k in seq_along(ev)) out$row[[paste0("ev", k)]] <- ev[k]
  out$truth <- list(B_true = r$truth$B, sigma2_true = r$truth$sigma2,
                    mu_true = r$truth$mu, phi_true = r$truth$phi)
  rm(r)
  gc(verbose = FALSE)
  out
}

found <- list_result_files(res_dirs)
expected <- expand.grid(J = J_grid, Q = Q_grid, seed = seeds, method = methods,
                        stringsAsFactors = FALSE)
if (only_expected) found <- found[key_of(found) %in% key_of(expected), ]
cat(sprintf("reading %d files with %d core(s)\n", nrow(found), n_cores))

res <- parallel::mclapply(seq_len(nrow(found)), function(i) {
  read_one(found$path[i], found$J[i], found$Q[i], found$seed[i], found$method[i])
}, mc.cores = n_cores, mc.preschedule = FALSE)

is_good <- vapply(res, function(x) is.list(x) && is.data.frame(x$row), logical(1))
for (i in which(!is_good)) {
  res[[i]] <- list(row = data.frame(J = found$J[i], Q = found$Q[i], seed = found$seed[i],
                                    method = found$method[i], file_exists = TRUE,
                                    read_ok = FALSE,
                                    read_error = "worker failed (possible out-of-memory)",
                                    stringsAsFactors = FALSE),
                   arrays = list(), truth = NULL)
}

present <- bind_fill(lapply(res, function(x) x$row))
present$expected <- key_of(present) %in% key_of(expected)
present$dup_key <- duplicated(key_of(present)) | duplicated(key_of(present), fromLast = TRUE)
n_present <- nrow(present)

miss <- expected[!(key_of(expected) %in% key_of(present)), , drop = FALSE]
if (nrow(miss) > 0) {
  miss$file_exists <- FALSE
  miss$expected <- TRUE
  miss$dup_key <- FALSE
  est <- bind_fill(list(present, miss))
} else {
  est <- present
}
n_miss <- nrow(est) - n_present

ord <- order(est$method, est$J, est$Q, est$seed)
est <- est[ord, ]
for (f in array_fields) {
  col <- lapply(res, function(x) x$arrays[[f]])
  est[[f]] <- c(col, vector("list", n_miss))[ord]
}
rownames(est) <- NULL

first_cols <- c("J", "Q", "seed", "method", "expected", "file_exists", "read_ok", "read_error",
                "dup_key", "K", "ok", "converged", "sigma2", "sigma2_last", "obj_value",
                "elapsed_sec", metric_names, "sigma2_bias_pct",
                grep("^ev[0-9]+$", names(est), value = TRUE), "B_dim_ok", "u_frac")
first_cols <- intersect(first_cols, names(est))
scalar_cols <- setdiff(names(est), c(array_fields))
est <- est[, c(first_cols, setdiff(scalar_cols, first_cols), array_fields)]

has_truth <- vapply(res, function(x) !is.null(x$truth), logical(1))
keys_all <- key3_of(do.call(rbind, lapply(res, function(x) x$row[c("J", "Q", "seed")])))
grp <- split(which(has_truth), keys_all[has_truth])
same_truth <- function(a, b) identical(a$B_true, b$B_true) && identical(a$sigma2_true, b$sigma2_true)
truth_conflict <- vapply(grp, function(ii) {
  any(!vapply(ii[-1], function(j) same_truth(res[[ii[1]]]$truth, res[[j]]$truth), logical(1)))
}, logical(1))
first_idx <- vapply(grp, function(ii) ii[1], integer(1))
truth <- do.call(rbind, lapply(first_idx, function(i) res[[i]]$row[c("J", "Q", "seed")]))
truth$K <- vapply(first_idx, function(i) null2na(res[[i]]$row$K), numeric(1))
truth$sigma2_true <- vapply(first_idx, function(i) null2na(res[[i]]$truth$sigma2_true), numeric(1))
for (f in c("B_true", "mu_true", "phi_true")) {
  truth[[f]] <- lapply(first_idx, function(i) res[[i]]$truth[[f]])
}
ev_true <- lapply(first_idx, function(i) ev_of(res[[i]]$truth$B_true))
for (k in seq_len(max(lengths(ev_true)))) {
  truth[[paste0("ev", k)]] <- vapply(ev_true, function(v) if (length(v) >= k) v[k] else NA_real_, numeric(1))
}
truth <- truth[order(truth$J, truth$Q, truth$seed), ]
rownames(truth) <- NULL

# summary(truth)

out_file <- file.path(out_dir, "estimates_coap.rds")
saveRDS(list(estimates = est, truth = truth), out_file)
write.csv(est[, !(names(est) %in% array_fields)], file.path(out_dir, "estimates_coap.csv"), row.names = FALSE)

cat("\nfile status by method and Q\n")
nkeys <- function(d, cond) length(unique(key_of(d[cond, , drop = FALSE])))
status_tbl <- do.call(rbind, lapply(
  split(est, list(est$method, est$Q), drop = TRUE),
  function(d) data.frame(method = d$method[1], Q = d$Q[1],
                         n_expected = nkeys(d, d$expected),
                         n_found = nkeys(d, d$expected & d$file_exists),
                         n_read_ok = nkeys(d, d$expected & d$read_ok %in% TRUE),
                         n_stray = sum(!d$expected))))
rownames(status_tbl) <- NULL
print(status_tbl)

cat("\nelapsed_sec by method and Q\n")
ok_rows <- est[est$read_ok %in% TRUE, ]
fstat <- function(x, f) if (all(is.na(x))) NA else round(f(x, na.rm = TRUE), 2)
time_tbl <- do.call(rbind, lapply(
  split(ok_rows, list(ok_rows$method, ok_rows$Q), drop = TRUE),
  function(d) data.frame(method = d$method[1], Q = d$Q[1], n = nrow(d),
                         n_time = sum(!is.na(d$elapsed_sec)),
                         min = fstat(d$elapsed_sec, min),
                         median = fstat(d$elapsed_sec, median),
                         mean = fstat(d$elapsed_sec, mean),
                         max = fstat(d$elapsed_sec, max))))
rownames(time_tbl) <- NULL
print(time_tbl)

cat("\nchecks\n")
cat("read failures:", sum(est$read_ok %in% FALSE), "\n")
cat("read ok but ok = FALSE:", sum(ok_rows$ok %in% FALSE), "\n")
cat("B dimension mismatches:", sum(ok_rows$B_dim_ok %in% FALSE), "\n")
cat("duplicated keys:", nkeys(est, est$dup_key %in% TRUE), "\n")
cat("stray files (not in grid):", sum(!est$expected), "\n")
cat("(J,Q,seed) with conflicting truth:", sum(truth_conflict), "\n")
if (any(truth_conflict)) print(names(truth_conflict)[truth_conflict])
cat("max u_frac by method (gauge, should be near 0)\n")
print(tapply(ok_rows$u_frac, ok_rows$method, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE)))
cat(sprintf("\nestimates rows %d, truth rows %d, saved %.1f MB\n",
            nrow(est), nrow(truth), file.info(out_file)$size / 2^20))
