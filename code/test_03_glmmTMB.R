library(glmmTMB) # for rr() comparison on the level of Poisson surrogate
library(bench) # timing
library(kableExtra)

stopifnot(packageVersion("glmmTMB") >= "1.1.8")   # covariance rr() requires

source("pfa_woodbury.R")
source("sim_data.R")
source("basic_functions.R")

# Check idiosyncratic term
nlevels(long$obs) == length(unique(interaction(long$group, long$category, drop = TRUE)))

# fit glmmTMB reduced-rank reference
fit_glmmTMB_ref <- function(Y, X, group, K, time_limit = Inf) {
  long <- build_long(Y, X, group)
  form <- as.formula(paste0("count ~ category + rr(category + 0 | group, d = ", K, ") + (1 | obs)")) 
  t0 <- proc.time()["elapsed"]
  fit <- tryCatch(
    glmmTMB(form, data = long, family = poisson(link = "log"),
            offset = long$log_total,
            control = glmmTMBControl(
              optCtrl = list(iter.max = 1000, eval.max = 1000))),
    error = function(e) NULL)
  el <- proc.time()["elapsed"] - t0
  if (is.null(fit)) return(list(B = NULL, time = el, ok = FALSE))
  
  conv <- isTRUE(fit$sdr$pdHess)
  
  # extract rank-K loadings from the rr covariance (compare the K eigen only because of B's indeterminacy)
  vc <- VarCorr(fit)$cond$group
  e <- eigen((vc + t(vc)) / 2, symmetric = TRUE)
  L <- e$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(e$values[1:K], 0)), K)
  s2_gt <- tryCatch(sigma(fit)^2, error = function(e) NA)
  list(B = L, time = el, ok = TRUE)
}

cat(sprintf("%5s %5s %5s | %14s %14s %16s\n",
            "Q", "K", "J", "wb vs true", "glmmTMB vs true", "wb vs glmmTMB"))

for (Q in c(20, 30, 50)) {
  K <- 2; J <- 30; Nj <- 15
  dat <- simulate_pfa_data(Q = Q, K = K, J = J, N_per_group = Nj, seed = 11)
  
  fw <- fit_pfa_woodbury(dat$Y, dat$X, dat$group, K = K, M = dat$M,
                         max_iter = 30, verbose = FALSE, exact_Shat = TRUE)
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

cat(sprintf("%6s | %14s %14s %14s %14s %10s %10s %8s\n",
            "Q", "Woodbury (s)", "wb mem (MB)", "glmmTMB (s)", "gT mem (MB)",
            "wb vs gT", "gt vs tru", "conv?"))
Qs <- c(20, 50, 100, 150)
results <- data.frame(
  Q = numeric(),
  Woodbury_Time = numeric(),
  Woodbury_Memory = numeric(),
  glmmTMB_Time = numeric(),
  glmmTMB_Memory = numeric(),
  Subspace_Distance = numeric(),
  glmmTMB_vs_True = numeric(),
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
  t_wb  <- as.numeric(b_wb$median)
  m_wb  <- as.numeric(b_wb$mem_alloc) / 1024^2
  
  b_gt <- tryCatch(
    bench::mark(fr <<- fit_glmmTMB_ref(dat$Y, dat$X, dat$group, K = K),
                iterations = 2, check = FALSE, memory = TRUE), 
    error = function(e) NULL)
  ok <- !is.null(b_gt) && fr$ok
  t_gt <- if (ok) as.numeric(b_gt$median) else NA
  m_gt <- if (ok) as.numeric(b_gt$mem_alloc) / 1024^2 else NA
  d <- if (ok) subspace_dist(fw$B, fr$B) else NA
  d_gt <- if (ok) subspace_dist(fr$B, dat$true$B) else NA
  conv <- if (ok) isTRUE(fr$converged) else FALSE
  
  cat(sprintf("%6d | %14.2f %14.1f %14s %14s %10s %10s %8s\n",
              Q, t_wb, m_wb,
              if (ok) sprintf("%.2f", t_gt) else "FAILED/slow",
              if (ok) sprintf("%.1f", m_gt) else "-",
              if (is.na(d)) "-" else sprintf("%.4f", d),
              if (is.na(d_gt)) "-" else sprintf("%.4f", d_gt),
              if (ok) ifelse(conv, "yes", "NO") else "-"))
  results <- rbind(
    results,
    data.frame(
      Q = Q,
      "Woodbury Time (s)" = t_wb,
      "Woodbury_Memory (MB)" = m_wb,
      "glmmTMB_Time (S)" = t_gt,
      "glmmTMB_Memory (MB)" = m_gt,
      "Subspace_Distance" = d,
      "glmmTMB_vs_True" = d_gt,
      "Converged" = ifelse(ok, ifelse(conv, "yes", "no"), "failed")
    )
  )
}
latex_table <- kbl(
  results,
  format = "latex",
  digits = 2,
  booktabs = TRUE,
  caption = "Woodbury vs glmmTMB.",
  label = "tab:benchmark"
) |>
  kable_styling(latex_options = "hold_position")
writeLines(latex_table, "wd_glmm_conmparison.tex")