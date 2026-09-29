suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

default_key_metrics <- c("B_dist", "tucker", "ang_max", "rv", "rmse", "sig_err")
default_logev_candidates <- c("log_evidence", "log_ev", "log_evidence_trace", "log_ev_trace")

is_scal <- function(z) {
  (is.numeric(z) || is.logical(z) || inherits(z, "difftime")) &&
    length(z) == 1 && is.null(dim(z))
}

as_num1 <- function(z) {
  if (is.null(z)) return(NA_real_)
  if (inherits(z, "difftime")) return(as.numeric(z, units = "secs"))
  if (is_scal(z)) return(as.numeric(z))
  NA_real_
}

cfg_str <- function(x) if (is.null(x)) NA_character_ else paste(x, collapse = ",")

make_tag <- function(J, Q, method) {
  sprintf("J%d_Q%d_%s", as.integer(J), as.integer(Q), method)
}

seed_file_pattern <- function(J, Q, method) {
  sprintf("^J%d_Q%d_seed([0-9]+)_%s\\.rds$", as.integer(J), as.integer(Q), method)
}

list_seed_files <- function(res_dir, J, Q, method) {
  pattern <- seed_file_pattern(J, Q, method)
  files <- list.files(res_dir, pattern = pattern, full.names = TRUE)
  seeds <- as.integer(sub(pattern, "\\1", basename(files)))
  ord <- order(seeds)
  list(files = files[ord], seeds = seeds[ord], pattern = pattern)
}

describe_fit <- function(path) {
  r1 <- readRDS(path)
  str(r1$config)
  print(names(r1$truth))
  fit_info <- tibble(
    field = names(r1$fit),
    class = map_chr(r1$fit, ~ class(.x)[1]),
    length = map_int(r1$fit, length),
    dim = map_chr(r1$fit, ~ paste(dim(.x), collapse = "x")),
    MB = round(map_dbl(r1$fit, ~ as.numeric(object.size(.x)) / 2^20), 3)
  )
  print(fit_info, n = Inf)
  invisible(fit_info)
}

extract_one <- function(path, seed_fn, logev_candidates = default_logev_candidates) {
  r <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(r, "error")) {
    return(list(row = tibble(seed = seed_fn, read_ok = FALSE,
                             read_msg = conditionMessage(r)),
                traces = NULL))
  }
  cf <- r$config
  tr <- r$truth
  fit <- r$fit

  lev_name <- intersect(logev_candidates, names(fit))[1]
  tr_names <- grep("_trace$", names(fit), value = TRUE)
  if (!is.na(lev_name)) tr_names <- union(lev_name, tr_names)
  tr_names <- tr_names[map_lgl(fit[tr_names], ~ (is.numeric(.x) || is.logical(.x)) &&
                                 is.null(dim(.x)) && length(.x) > 0)]

  traces <- map_dfr(tr_names, function(nm) {
    v <- as.numeric(fit[[nm]])
    tibble(metric = sub("_trace$", "", nm), iter = seq_along(v), value = v)
  })

  lev <- if (!is.na(lev_name)) as.numeric(fit[[lev_name]]) else numeric(0)
  peak_iter <- if (length(lev) > 0 && any(is.finite(lev))) which.max(lev) else NA_integer_
  n_drop <- if (length(lev) > 1) {
    d <- diff(lev)
    sum(d < -1e-8 * pmax(1, abs(head(lev, -1))), na.rm = TRUE)
  } else NA_integer_
  max_drop <- if (length(lev) > 1) max(c(0, -diff(lev)), na.rm = TRUE) else NA_real_

  fin_wide <- traces %>%
    group_by(metric) %>%
    summarise(final = value[n()],
              peak = if (!is.na(peak_iter) && n() >= peak_iter) value[peak_iter] else NA_real_,
              .groups = "drop") %>%
    pivot_wider(names_from = metric, values_from = c(final, peak),
                names_glue = "{metric}_{.value}")

  sc <- map_dbl(fit[map_lgl(fit, is_scal)], as_num1)

  row <- tibble(
    seed = seed_fn,
    read_ok = TRUE,
    read_msg = NA_character_,
    cfg_J = cf$J %||% NA,
    cfg_Q = cf$Q %||% NA,
    cfg_K = cf$K %||% NA,
    cfg_P = cf$P %||% NA,
    cfg_seed = cf$seed %||% NA,
    cfg_method = cf$method %||% NA_character_,
    cfg_M_range = cfg_str(cf$M_range),
    cfg_N_range = cfg_str(cf$N_per_group_range),
    cfg_sigma2_in = cf$sigma2_true_input %||% NA,
    truth_fp = signif(sum(tr$B^2) + sum(tr$mu * seq_along(tr$mu)), 12),
    sigma2_true = as_num1(tr$sigma2),
    sigma2_hat = as_num1(fit$sigma2),
    B_finite = !is.null(fit$B) && all(is.finite(fit$B)),
    n_iter_fit = as_num1(fit$iterations),
    logev_field = lev_name,
    peak_iter = peak_iter,
    n_logev_drop = n_drop,
    max_logev_drop = max_drop,
    stop_reason = if (is.character(fit$stop_reason)) fit$stop_reason[1] else NA_character_,
    total_time_sec = if (is.null(fit$iter_wall_time_trace)) NA_real_ else sum(fit$iter_wall_time_trace, na.rm = TRUE)
  )
  if (length(sc) > 0) {
    row <- bind_cols(row, as_tibble(as.list(sc)) %>% rename_with(~ paste0("fit_", .x)))
  }
  if (ncol(fin_wide) > 0) row <- bind_cols(row, fin_wide)

  list(row = row, traces = traces %>% mutate(seed = seed_fn, .before = 1))
}

mode_chr <- function(x) names(sort(table(as.character(x)), decreasing = TRUE))[1]

same_as_mode <- function(x) as.character(x) %in% mode_chr(x)

rob_z <- function(x, floor_rel = 1e-6) {
  m <- median(x, na.rm = TRUE)
  s <- mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- sd(x, na.rm = TRUE)
  s <- max(c(s, floor_rel * max(1, abs(m))), na.rm = TRUE)
  (x - m) / s
}

tr_stat <- function(hH, H) {
  M <- t(H) %*% hH %*% MASS::ginv(t(hH) %*% hH) %*% t(hH) %*% H
  sum(diag(M)) / sum(diag(t(H) %*% H))
}

cc <- function(M) sweep(M, 2, colMeans(M))

centering_matrix <- function(Q) diag(Q) - matrix(1 / Q, Q, Q)

oracle_eigen <- function(lam) {
  Q <- nrow(lam)
  J <- ncol(lam)
  P1 <- centering_matrix(Q)
  S <- P1 %*% lam %*% t(lam) %*% P1 / J
  list(ev = eigen(S, symmetric = TRUE), P1 = P1, Q = Q, J = J)
}

oracle_fa <- function(ev, K, Q) {
  s2 <- mean(ev$values[(K + 1):(Q - 1)])
  B <- ev$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(ev$values[1:K] - s2, 0)), K)
  list(B = B, sigma2 = s2)
}

align_lambda_hat <- function(lh, lam) {
  if (cor(c(t(lh)), c(lam)) > cor(c(lh), c(lam))) lh <- t(lh)
  lh
}

factor_scores <- function(B, s2, P1, lam) {
  K <- ncol(B)
  t(solve(crossprod(B) + s2 * diag(K), crossprod(B, P1 %*% lam)))
}

posthoc_one <- function(path, seed) {
  r <- readRDS(path)
  fit <- r$fit
  lam <- r$truth$lambda
  Bt <- r$truth$B
  Ft <- r$truth$F
  oe <- oracle_eigen(lam)
  P1 <- oe$P1
  Q <- oe$Q
  J <- oe$J

  or_cfg <- oracle_fa(oe$ev, r$config$K, Q)
  Bo_plt <- apply_PLT(or_cfg$B)

  or_B <- oracle_fa(oe$ev, ncol(Bt), Q)
  Fo <- factor_scores(or_B$B, or_B$sigma2, P1, lam)

  lh <- align_lambda_hat(fit$lambda_hat, lam)
  Fh <- factor_scores(fit$B, fit$sigma2, P1, lh)

  tibble(
    seed = seed,
    oracle_B_dist = subspace_distance(Bo_plt, Bt),
    oracle_tucker = tucker_congruence(Bo_plt, Bt)$mean,
    oracle_ang_max = max(principal_angles(Bo_plt, Bt)),
    oracle_sigma2 = or_cfg$sigma2,
    Tr_B_wb = tr_stat(fit$B, Bt),
    Tr_B_oracle = tr_stat(or_B$B, Bt),
    Tr_H_wb = tr_stat(Fh, Ft),
    Tr_H_oracle = tr_stat(Fo, Ft),
    lh_group_mean_share = J * sum(rowMeans(lh)^2) / sum(lh^2),
    Fh_mean_share = J * sum(colMeans(Fh)^2) / sum(Fh^2),
    loss_center = J * sum(colMeans(Ft)^2) / sum(Ft^2),
    Tr_H_wb_c = tr_stat(cc(Fh), cc(Ft)),
    Tr_H_oracle_c = tr_stat(cc(Fo), cc(Ft))
  )
}

check_fadmr_funs <- function(funs = c("apply_PLT", "subspace_distance",
                                      "tucker_congruence", "principal_angles")) {
  miss <- funs[!vapply(funs, function(f) exists(f, mode = "function"), logical(1))]
  if (length(miss) > 0) {
    stop("load the FA-DMR helper functions first, missing: ", paste(miss, collapse = ", "))
  }
  invisible(TRUE)
}

collect_grid <- function(out_dir, method, prefix, element) {
  pat <- sprintf("^%s_J[0-9]+_Q[0-9]+_%s\\.rds$", prefix, method)
  fs <- list.files(out_dir, pattern = pat, full.names = TRUE)
  bind_rows(map(fs, ~ readRDS(.x)[[element]]))
}
