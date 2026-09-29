res_dirs <- c("G:/J100_Q_seed10_COAP")
out_dir <- c("G:/J100_Q_seed10_COAP_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

J_grid <- c(100)
Q_grid <- c(50, 100, 150, 200, 300, 400)
seeds <- 1:100
methods <- c("coap_obs", "coap_group")
n_cores <- max(1L, as.integer(Sys.getenv("PBS_NCPUS", "1")))

null2na <- function(x) if (is.null(x) || length(x) != 1) NA else x

coap_field <- function(fit, name, alt = name) {
  if (!is.null(fit$best)) fit$best[[name]] else if (!is.null(fit[[name]])) fit[[name]] else fit[[alt]]
}

fit_views <- list(
  wb = function(fit) list(
    ok = !is.null(fit$B) && all(is.finite(fit$B)),
    converged = isTRUE(fit$converged),
    iterations = fit$iterations,
    stop_reason = fit$stop_reason,
    best_iter = fit$best_iter,
    elapsed_sec = fit$meta$elapsed_sec,
    pbs_job_id = fit$meta$pbs_job_id,
    has_trace = !is.null(fit$B_dist_trace),
    error = NA_character_),
  oracle = function(fit) list(
    ok = fit$ok,
    converged = fit$converged,
    elapsed_sec = fit$meta$elapsed_sec,
    error = fit$error),
  glmmTMB = function(fit) list(
    ok = fit$best$ok,
    converged = fit$best$converged,
    elapsed_sec = fit$meta$elapsed_sec,
    pbs_job_id = fit$meta$pbs_job_id,
    error = fit$best$error,
    n_restart = length(fit$all_restarts),
    best_idx = fit$best_idx,
    best_restart_sec = fit$timing$best_restart_sec,
    all_restarts_sum_sec = fit$timing$all_restarts_sum_sec),
  coap_group = function(fit) list(
    ok = isTRUE(coap_field(fit, "ok")),
    converged = isTRUE(coap_field(fit, "converged")),
    hit_cap = coap_field(fit, "hit_cap"),
    elapsed_sec = fit$meta$elapsed_sec,
    pbs_job_id = fit$meta$pbs_job_id,
    error = fit$error
  ),
  coap_obs = function(fit) list(
    ok = isTRUE(coap_field(fit, "ok")),
    converged = isTRUE(coap_field(fit, "converged")),
    n_iter = coap_field(fit, "n_iter", alt = "iterations"),
    elapsed_sec = fit$meta$elapsed_sec,
    pbs_job_id = fit$meta$pbs_job_id,
    error = fit$error,
  )
)

default_view <- function(fit) list(
  ok = fit$ok,
  converged = fit$converged,
  iterations = fit$iterations,
  stop_reason = fit$stop_reason,
  elapsed_sec = fit$meta$elapsed_sec,
  error = fit$error)

truth_fp <- function(B) {
  if (is.null(B)) return(NA_character_)
  sprintf("%.10e_%.10e", sum(B^2), sum(B * (seq_along(B) %% 13)))
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
  info <- file.info(path)
  row <- data.frame(J = J, Q = Q, seed = seed, method = method,
                    path = normalizePath(path), file_exists = TRUE,
                    file_size_mb = round(info$size / 2^20, 3),
                    mtime = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
                    stringsAsFactors = FALSE)
  r <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(r, "error")) {
    row$read_ok <- FALSE
    row$read_error <- conditionMessage(r)
    return(row)
  }
  row$read_ok <- TRUE
  row$read_error <- NA_character_
  cfg <- r$config
  tr <- r$truth
  row$key_match <- isTRUE(cfg$J == J) && isTRUE(cfg$Q == Q) &&
    isTRUE(cfg$seed == seed) && isTRUE(cfg$method == method)
  row$K <- null2na(cfg$K)
  row$P <- null2na(cfg$P)
  row$N_lo <- null2na(cfg$N_per_group_range[1])
  row$N_hi <- null2na(cfg$N_per_group_range[2])
  row$M_lo <- null2na(cfg$M_range[1])
  row$M_hi <- null2na(cfg$M_range[2])
  row$sigma2_true_input <- null2na(cfg$sigma2_true_input)
  row$N_total <- if (is.null(tr$n_per_group)) NA else sum(tr$n_per_group)
  row$n_group_min <- if (is.null(tr$n_per_group)) NA else min(tr$n_per_group)
  row$n_group_max <- if (is.null(tr$n_per_group)) NA else max(tr$n_per_group)
  row$sigma2_true <- null2na(tr$sigma2)
  row$truth_fp <- truth_fp(tr$B)
  view_fn <- if (is.null(fit_views[[method]])) default_view else fit_views[[method]]
  view <- tryCatch(view_fn(r$fit),
                   error = function(e) list(error = paste("view failed:", conditionMessage(e))))
  for (nm in names(view)) row[[nm]] <- null2na(view[[nm]])
  rm(r)
  gc(verbose = FALSE)
  row
}

found <- list_result_files(res_dirs)
# expected <- expand.grid(J = J_grid, Q = Q_grid, seed = seeds, method = methods,
#                         stringsAsFactors = FALSE)
# found <- found[key_of(found) %in% key_of(expected), ]
cat(sprintf("reading %d files with %d core(s)\n", nrow(found), n_cores))

rows <- parallel::mclapply(seq_len(nrow(found)), function(i) {
  read_one(found$path[i], found$J[i], found$Q[i], found$seed[i], found$method[i])
}, mc.cores = n_cores, mc.preschedule = FALSE)

failed <- which(!vapply(rows, is.data.frame, logical(1)))
for (i in failed) {
  rows[[i]] <- data.frame(J = found$J[i], Q = found$Q[i], seed = found$seed[i],
                          method = found$method[i], path = normalizePath(found$path[i]),
                          file_exists = TRUE, read_ok = FALSE,
                          read_error = "worker failed (possible out-of-memory)",
                          stringsAsFactors = FALSE)
}

present <- bind_fill(rows)
expected <- expand.grid(J = J_grid, Q = Q_grid, seed = seeds, method = methods,
                        stringsAsFactors = FALSE)
present$expected <- key_of(present) %in% key_of(expected)
present$dup_key <- duplicated(key_of(present)) | duplicated(key_of(present), fromLast = TRUE)

miss <- expected[!(key_of(expected) %in% key_of(present)), , drop = FALSE]
if (nrow(miss) > 0) {
  miss$path <- NA_character_
  miss$file_exists <- FALSE
  miss$expected <- TRUE
  miss$dup_key <- FALSE
  manifest <- bind_fill(list(present, miss))
} else {
  manifest <- present
}

first_cols <- c("J", "Q", "seed", "method", "expected", "file_exists", "read_ok", "read_error",
                "dup_key", "key_match", "path", "file_size_mb", "mtime",
                "K", "P", "N_lo", "N_hi", "M_lo", "M_hi", "sigma2_true_input",
                "N_total", "n_group_min", "n_group_max", "sigma2_true", "truth_fp",
                "ok", "converged", "iterations", "stop_reason", "elapsed_sec", "pbs_job_id", "error")
first_cols <- intersect(first_cols, names(manifest))
manifest <- manifest[, c(first_cols, setdiff(names(manifest), first_cols))]
manifest <- manifest[order(manifest$method, manifest$J, manifest$Q, manifest$seed), ]
rownames(manifest) <- NULL
manifest[c("path", "mtime", "pbs_job_id", "sigma2_true_input")] <- NULL

summary(manifest)

saveRDS(manifest, file.path(out_dir, "manifest_coap.rds"))
write.csv(manifest, file.path(out_dir, "manifest_coap.csv"), row.names = FALSE)

cat("\nfile status by method and Q\n")
nkeys <- function(d, cond) length(unique(key_of(d[cond, , drop = FALSE])))
status_tbl <- do.call(rbind, lapply(
  split(manifest, list(manifest$method, manifest$Q), drop = TRUE),
  function(d) data.frame(method = d$method[1], Q = d$Q[1],
                         n_expected = nkeys(d, d$expected),
                         n_found = nkeys(d, d$expected & d$file_exists),
                         n_read_ok = nkeys(d, d$expected & d$read_ok %in% TRUE),
                         n_stray = sum(!d$expected),
                         size_gb = round(sum(d$file_size_mb, na.rm = TRUE) / 1024, 3))))
rownames(status_tbl) <- NULL
print(status_tbl)

ok_rows <- manifest[manifest$read_ok %in% TRUE, ]

cat("\nconfiguration groups (each method should show one group)\n")
cfg_key <- with(ok_rows, paste0("K=", K, " P=", P, " N=[", N_lo, ",", N_hi, "] M=[", M_lo, ",", M_hi,
                                "] sigma2_in=", sigma2_true))
cfg_tbl <- as.data.frame(table(method = ok_rows$method, config = cfg_key), stringsAsFactors = FALSE)
print(cfg_tbl[cfg_tbl$Freq > 0, ], row.names = FALSE)

cat("\nstop_reason counts by method and Q\n")
sr_tbl <- as.data.frame(table(method = ok_rows$method, Q = ok_rows$Q,
                              stop_reason = ok_rows$stop_reason, useNA = "ifany"),
                        stringsAsFactors = FALSE)
print(sr_tbl[sr_tbl$Freq > 0, ], row.names = FALSE)

cat("\niterations by method and Q\n")
fstat <- function(x, f) if (all(is.na(x))) NA else f(x, na.rm = TRUE)
iter_tbl <- do.call(rbind, lapply(
  split(ok_rows, list(ok_rows$method, ok_rows$Q), drop = TRUE),
  function(d) data.frame(method = d$method[1], Q = d$Q[1], n = nrow(d),
                         iter_min = fstat(d$iterations, min),
                         iter_median = fstat(d$iterations, median),
                         iter_max = fstat(d$iterations, max))))
rownames(iter_tbl) <- NULL
print(iter_tbl)

cat("\nchecks\n")
cat("read failures:", sum(manifest$read_ok %in% FALSE), "\n")
cat("filename vs config mismatches:", sum(ok_rows$key_match %in% FALSE), "\n")
cat("duplicated keys:", nkeys(manifest, manifest$dup_key %in% TRUE), "\n")
cat("stray files (not in grid):", sum(!manifest$expected), "\n")
fp <- ok_rows[!is.na(ok_rows$truth_fp), ]
n_fp <- tapply(paste(fp$truth_fp, fp$N_total), paste(fp$J, fp$Q, fp$seed), function(x) length(unique(x)))
cat("keys with differing truth:", sum(n_fp > 1), "\n")
if (any(n_fp > 1)) print(names(n_fp)[n_fp > 1])
