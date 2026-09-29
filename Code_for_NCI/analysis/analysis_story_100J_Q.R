source("analysis_utils.R")

if (!exists("rr")) rr <- readRDS("estimates.rds")
if (!exists("mf")) mf <- readRDS("manifest.rds")

fig_dir <- "figs_story"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

step_title <- function(txt) {
  cat("\n", strrep("=", 72), "\n", txt, "\n", strrep("=", 72), "\n", sep = "")
}


as_est_df <- function(x, nm) {
  need <- c("J", "Q", "seed", "method")
  if (is.data.frame(x)) return(as_tibble(x))
  if (is.list(x) && all(need %in% names(x)) && length(unique(lengths(x[need]))) == 1) {
    return(as_tibble(x[lengths(x) == length(x[[need[1]]])]))
  }
  hit <- keep(x, ~ is.data.frame(.x) && all(need %in% names(.x)))
  if (length(hit) >= 1) return(as_tibble(bind_rows(hit)))
  rows <- keep(x, ~ is.list(.x) && all(need %in% names(.x)))
  if (length(rows) >= 1) return(bind_rows(map(rows, ~ as_tibble(keep(.x, ~ is.atomic(.x) && length(.x) == 1)))))
  stop(nm, " is a ", class(x)[1], " with elements: ", paste(head(names(x), 20), collapse = ", "))
}

key <- c("J", "Q", "seed", "method")
rr_df <- as_est_df(rr, "rr")
mf_df <- as_est_df(mf, "mf")
cat(sprintf("rr coerced to %d x %d, mf coerced to %d x %d\n",
            nrow(rr_df), ncol(rr_df), nrow(mf_df), ncol(mf_df)))
mf2 <- mf_df %>%
  mutate(trB2_true = as.numeric(sub("_.*$", "", truth_fp))) %>%
  select(all_of(key), sigma2_true, truth_fp, trB2_true, N_total, n_group_min, n_group_max,
         N_lo, N_hi, M_lo, M_hi, P, expected, file_exists, read_ok, key_match, dup_key,
         iterations, stop_reason, best_iter, has_trace)

rrj <- rr_df %>%
  select(-any_of(setdiff(names(mf2), key))) %>%
  inner_join(mf2, by = key)

J0 <- unique(rrj$J)
K0 <- unique(rrj$K)
stopifnot(length(J0) == 1, length(K0) == 1)

save_fig <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(fig_dir, sprintf("J%d_%s.png", as.integer(J0), name)), p,
         width = w, height = h, dpi = 150)
  print(p)
  invisible(p)
}

theme_set(theme_bw(base_size = 12))
col_m <- c(oracle = "grey55", wb = "#2c7fb8")
mcse <- function(x) sd(x) / sqrt(length(x))

spike_from_ev <- function(ev, s2, g) {
  l <- ev + s2
  a <- l - s2 * (1 + g)
  (a + sqrt(pmax(a^2 - 4 * g * s2^2, 0))) / 2
}

cos2_spiked <- function(th, s2, g) (1 - g * s2^2 / th^2) / (1 + g * s2 / th)

step_title("STEP 0. Data integrity")

cat(sprintf("rows in manifest %d, rows after join with estimates %d\n", nrow(mf), nrow(rrj)))
stopifnot(nrow(rrj) == nrow(mf))
flag_tab <- rrj %>%
  summarise(expected = all(expected), file_exists = all(file_exists), read_ok = all(read_ok),
            key_match = all(key_match), no_dup = !any(dup_key))
print(flag_tab)
stopifnot(all(unlist(flag_tab)))

fp_chk <- rrj %>% group_by(Q, seed) %>% summarise(n_fp = n_distinct(truth_fp), .groups = "drop")
stopifnot(all(fp_chk$n_fp == 1))
cat("truth_fp is identical for oracle and wb within every (Q, seed). The pairing is on the same truth.\n")

chk <- rrj %>%
  group_by(Q, method) %>%
  summarise(n = n(), n_seed = n_distinct(seed), all_ok = all(ok),
            all_B_dim_ok = all(B_dim_ok), max_u_frac = max(u_frac),
            n_na_d = sum(is.na(d_true)), n_na_s2 = sum(is.na(sigma2)), .groups = "drop")
print(chk, n = Inf)

pair_chk <- rrj %>% count(Q, seed)
stopifnot(all(pair_chk$n == 2), all(chk$all_ok), all(chk$all_B_dim_ok))
cat("Every (Q, seed) has exactly one oracle and one wb fit.\n")
cat(sprintf("Largest energy in the 1_Q gauge direction: %.2e (should be ~0).\n", max(rrj$u_frac)))


step_title("STEP 1. Design")

s2_true <- unique(rrj$sigma2_true)
stopifnot(length(s2_true) == 1)
Q_vals <- sort(unique(rrj$Q))
cat(sprintf("J = %d, K = %d, P = %d, sigma2_true = %.4f\n",
            as.integer(J0), as.integer(K0), as.integer(unique(rrj$P)), s2_true))
cat(sprintf("N per group in [%d, %d], M per observation in [%d, %d]\n",
            as.integer(unique(rrj$N_lo)), as.integer(unique(rrj$N_hi)),
            as.integer(unique(rrj$M_lo)), as.integer(unique(rrj$M_hi))))
cat("Q values:", Q_vals, "\n")
cat("Seeds per (Q, method):", unique(chk$n_seed), "\n")
print(rrj %>% filter(method == "oracle") %>% group_by(Q) %>%
        summarise(N_total_mean = mean(N_total), N_total_sd = sd(N_total),
                  trB2_true_mean = mean(trB2_true), .groups = "drop"), n = Inf)

rr2 <- rrj %>% mutate(Qf = factor(Q, levels = Q_vals),
                     method = factor(method, levels = c("oracle", "wb")))

w <- rr2 %>%
  select(Q, Qf, seed, method, trB2_true, N_total, sigma2, sigma2_bias_pct, d_true, ang_max, ev1, ev2) %>%
  pivot_wider(names_from = method,
              values_from = c(sigma2, sigma2_bias_pct, d_true, ang_max, ev1, ev2))


step_title("STEP 2. truth -> oracle, sigma2")

cat("The oracle sees the true lambda_j. It is PPCA on J vectors in Q-1 identified dimensions.\n")
cat("Theory. The top K sample eigenvalues absorb about gamma*sigma2 each, gamma = (Q-1)/J.\n")
cat("So E[sigma2_oracle]/sigma2 - 1 = -K(Q-1)/(J(Q-1-K)), close to -K/J.\n")
cat("The residual eigenvalues average Q-1-K chi-square pieces, so sd = sigma2*sqrt(2/(J(Q-1-K))).\n\n")

or_s2 <- w %>%
  mutate(m_eff = Q - 1 - K0,
         pred_bias = -100 * K0 * (Q - 1) / (J0 * m_eff),
         pred_sd = 100 * sqrt(2 / (J0 * m_eff)))

tab_or_s2 <- or_s2 %>%
  group_by(Q) %>%
  summarise(obs_bias = mean(sigma2_bias_pct_oracle), mcse = mcse(sigma2_bias_pct_oracle),
            pred_bias = first(pred_bias), z = (obs_bias - pred_bias) / mcse,
            obs_sd = sd(sigma2_bias_pct_oracle), pred_sd = first(pred_sd),
            sd_ratio = obs_sd / pred_sd, .groups = "drop")
print(tab_or_s2 %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf)
cat(sprintf("Reference -K/J = %.2f%%. |z| < 2 means theory and simulation agree within Monte Carlo error.\n",
            -100 * K0 / J0))

band <- or_s2 %>% distinct(Qf, pred_bias, pred_sd)
p2 <- ggplot(or_s2, aes(Qf, sigma2_bias_pct_oracle)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(fill = col_m["oracle"], alpha = 0.5, outlier.size = 0.8) +
  geom_errorbar(data = band, aes(x = Qf, ymin = pred_bias - 2 * pred_sd, ymax = pred_bias + 2 * pred_sd),
                inherit.aes = FALSE, width = 0.25, colour = "firebrick") +
  geom_point(data = band, aes(x = Qf, y = pred_bias), inherit.aes = FALSE,
             colour = "firebrick", size = 2.5) +
  labs(x = "Q", y = "oracle sigma2 bias (%)",
       title = "truth -> oracle: sigma2",
       subtitle = "red = theory mean -K(Q-1)/(J(Q-1-K)) with +/- 2 theory sd")
save_fig(p2, "step2_oracle_sigma2")


step_title("STEP 3. truth -> oracle, B")

cat("The chain has one more link than sigma2. truth B -> realised signal B (F'F/J) B' -> oracle.\n")
cat("The oracle sees lambda_j = B f_j + e_j, so it can only recover the signal in the realised factors.\n")
cat("theta_k = eigenvalues of B'B. The manifest gives tr(B'B) = theta_1 + theta_2 through truth_fp.\n")
cat("theta_k^real is recovered per seed from the oracle eigenvalues by inverting the spiked bias.\n\n")

or_B <- w %>%
  mutate(g = (Q - 1) / J0,
         th1 = spike_from_ev(ev1_oracle, sigma2_oracle, g),
         th2 = spike_from_ev(ev2_oracle, sigma2_oracle, g),
         th_sum = th1 + th2,
         share1 = th1 / th_sum,
         real_over_true = th_sum / trB2_true,
         pred_sd_real = sqrt(2 / J0) * sqrt(th1^2 + th2^2) / th_sum,
         omega_true = trB2_true / Q,
         omega1 = th1 / Q, omega2 = th2 / Q)

step_title("STEP 3a. truth -> realised signal")

cat("If f_j ~ N(0, I_K), the realised signal trace over the true trace has mean 1\n")
cat("and relative sd sqrt(2/J) * sqrt(theta_1^2 + theta_2^2) / (theta_1 + theta_2).\n\n")
tab_sig <- or_B %>%
  group_by(Q) %>%
  summarise(omega_true = mean(omega_true), ratio_mean = mean(real_over_true),
            ratio_mcse = mcse(real_over_true), ratio_sd = sd(real_over_true),
            pred_sd = mean(pred_sd_real), cor_true_real = cor(trB2_true, th_sum), .groups = "drop")
print(tab_sig %>% mutate(across(where(is.double), ~ round(.x, 4))), n = Inf)
cat("A mean near 1 says the spike recovery is unbiased against the true B.\n")
cat("The sd matching the prediction says the seed-level scatter is factor sampling, not an estimation error.\n")

p3a <- ggplot(or_B, aes(Qf, omega_true)) +
  geom_boxplot(fill = "#fdae6b", outlier.size = 0.8) +
  labs(x = "Q", y = "omega_1 + omega_2 = tr(B'B) / Q",
       title = "True signal per category is stable across Q",
       subtitle = "tr(B'B) from truth_fp in the manifest")
print(p3a)
save_fig(p3a, "step3a_omega_true_by_Q")

sd_band <- or_B %>% group_by(Qf) %>% summarise(s = mean(pred_sd_real), .groups = "drop")
p3b <- ggplot(or_B, aes(Qf, real_over_true)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_boxplot(fill = "#fdae6b", alpha = 0.6, outlier.size = 0.8) +
  geom_errorbar(data = sd_band, aes(x = Qf, ymin = 1 - 2 * s, ymax = 1 + 2 * s),
                inherit.aes = FALSE, width = 0.25, colour = "firebrick") +
  labs(x = "Q", y = "realised signal / true signal",
       title = "truth -> realised signal",
       subtitle = "red = 1 +/- 2 sd predicted from factor sampling alone")
save_fig(p3b, "step3b_realised_over_true")
print(p3b)

step_title("STEP 3b. realised signal -> oracle B")

cat("Spiked covariance theory gives cos2_k = (1 - gamma*s2^2/theta_k^2) / (1 + gamma*s2/theta_k).\n")
cat("This is the formula in bbp_floor(), with gamma = (Q-1)/J for the gauge-projected space.\n")
cat("Two versions are compared.\n")
cat("  truth-based. theta_k = tr(B'B)_true * share_k, s2 = sigma2_true. Uses the truth, ignores F sampling.\n")
cat("  realised. theta_k recovered from the oracle, s2 = sigma2_oracle. Conditions on the realised F.\n")
cat("share_k is taken from the oracle because the manifest stores only the trace of B'B.\n\n")

d_from_theta <- function(t1, t2, s2, g) {
  c1 <- cos2_spiked(t1, s2, g)
  c2 <- cos2_spiked(t2, s2, g)
  list(d = sqrt(2 * ((1 - c1) + (1 - c2))), ang = acos(sqrt(pmin(c1, c2))) * 180 / pi)
}

or_B <- or_B %>%
  mutate(d_pred_true = d_from_theta(trB2_true * share1, trB2_true * (1 - share1), s2_true, g)$d,
         ang_pred_true = d_from_theta(trB2_true * share1, trB2_true * (1 - share1), s2_true, g)$ang,
         d_pred = d_from_theta(th1, th2, sigma2_oracle, g)$d,
         ang_pred = d_from_theta(th1, th2, sigma2_oracle, g)$ang,
         d_lead = sqrt(2 * g * s2_true * (1 / th1 + 1 / th2)),
         d_ratio_pred = d_true_oracle / d_pred,
         d_ratio_pred_true = d_true_oracle / d_pred_true)

tab_or_B <- or_B %>%
  group_by(Q) %>%
  summarise(cor_seed_real = cor(d_true_oracle, d_pred),
            cor_seed_true = cor(d_true_oracle, d_pred_true),
            obs_over_pred_real = mean(d_ratio_pred),
            obs_over_pred_true = mean(d_ratio_pred_true),
            d_obs = mean(d_true_oracle), d_pred = mean(d_pred), d_pred_true = mean(d_pred_true),
            ang_obs = mean(ang_max_oracle), ang_pred = mean(ang_pred), .groups = "drop")
print(tab_or_B %>% mutate(across(where(is.double), ~ round(.x, 4))), n = Inf)
cat("Both versions should match in the mean. The realised version should track seeds more tightly.\n")
cat("Key point. To leading order 1 - cos2_k = gamma*s2/theta_k = s2/(J*omega_k) * (Q-1)/Q.\n")
cat("omega_k = theta_k/Q stays flat in Q, so Q cancels. The oracle error is set by J and omega.\n")

pred_pts <- or_B %>%
  group_by(Qf) %>%
  summarise(realised = mean(d_pred), truth_based = mean(d_pred_true), .groups = "drop") %>%
  pivot_longer(-Qf, names_to = "prediction", values_to = "d")
p3c <- ggplot(or_B, aes(Qf, d_true_oracle)) +
  geom_boxplot(fill = col_m["oracle"], alpha = 0.5, outlier.size = 0.8) +
  geom_point(data = pred_pts, aes(Qf, d, colour = prediction), size = 2.5,
             position = position_dodge(0.3)) +
  geom_line(data = pred_pts, aes(Qf, d, colour = prediction, group = prediction)) +
  scale_colour_manual(values = c(realised = "firebrick", truth_based = "darkgreen")) +
  labs(x = "Q", y = "oracle subspace distance to truth", colour = NULL,
       title = "realised signal -> oracle: B", subtitle = "no fitted constant")
save_fig(p3c, "step3c_oracle_d_true")
print(p3c)

sc_long <- or_B %>%
  select(Qf, seed, d_true_oracle, realised = d_pred, truth_based = d_pred_true) %>%
  pivot_longer(c(realised, truth_based), names_to = "prediction", values_to = "d_pred")
p3d <- ggplot(sc_long, aes(d_pred, d_true_oracle, colour = Qf)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.5, size = 1) +
  facet_wrap(~ prediction) +
  scale_colour_viridis_d() +
  labs(x = "predicted d", y = "observed oracle d_true", colour = "Q",
       title = "Seed-level check. Oracle vs spiked-model prediction")
save_fig(p3d, "step3d_oracle_pred_scatter", w = 10, h = 5)
print(p3d)

rat_long <- or_B %>%
  select(Qf, seed, realised = d_ratio_pred, truth_based = d_ratio_pred_true) %>%
  pivot_longer(-c(Qf, seed), names_to = "prediction", values_to = "ratio")
p3e <- ggplot(rat_long, aes(Qf, ratio, fill = prediction)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c(realised = "#fc9272", truth_based = "#a1d99b")) +
  labs(x = "Q", y = "observed / predicted", fill = NULL,
       title = "Oracle d_true relative to theory")
save_fig(p3e, "step3e_oracle_ratio_to_theory")
print(p3e)

step_title("STEP 4. oracle -> wb, paired comparison")

cat("Both methods use the same simulated data for a given seed, so every comparison is paired.\n")
cat("Pairing removes seed-to-seed variation in omega and isolates the algorithm layer.\n\n")

pw <- or_B %>%
  mutate(gap_s2_pct = sigma2_bias_pct_wb - sigma2_bias_pct_oracle,
         ratio_d = d_true_wb / d_true_oracle,
         ex_d2 = d_true_wb^2 - d_true_oracle^2,
         share_alg = ex_d2 / d_true_wb^2,
         gap_ang = ang_max_wb - ang_max_oracle,
         omega_min = pmin(omega1, omega2))

tab_pair <- pw %>%
  group_by(Q) %>%
  summarise(gap_s2_mcse = mcse(gap_s2_pct), gap_s2_pct = mean(gap_s2_pct),
            frac_wb_s2_lower = mean(sigma2_wb < sigma2_oracle),
            p_wilcox_s2 = wilcox.test(sigma2_wb, sigma2_oracle, paired = TRUE)$p.value,
            ratio_d_sd = sd(ratio_d), frac_wb_d_worse = mean(ratio_d > 1),
            ratio_d = mean(ratio_d),
            p_wilcox_d = wilcox.test(d_true_wb, d_true_oracle, paired = TRUE)$p.value,
            gap_ang = mean(gap_ang), .groups = "drop")
print(tab_pair %>% mutate(across(where(is.double), ~ signif(.x, 3))), n = Inf)

p4a <- ggplot(rr2, aes(Qf, sigma2_bias_pct, fill = method)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = -100 * K0 / J0, linetype = 3, colour = "firebrick") +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = col_m) +
  labs(x = "Q", y = "sigma2 bias (%)", fill = NULL,
       title = "sigma2: oracle stays at -K/J, wb drifts down with Q",
       subtitle = "dotted red line = -K/J")
save_fig(p4a, "step4a_sigma2_both")
print(p4a)

p4b <- ggplot(rr2, aes(Qf, d_true, fill = method)) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = col_m) +
  labs(x = "Q", y = "subspace distance to truth", fill = NULL,
       title = "B: oracle is flat in Q, wb grows with Q")
save_fig(p4b, "step4b_d_true_both")
print(p4b)

p4c <- ggplot(pw, aes(Qf, gap_s2_pct)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(fill = col_m["wb"], alpha = 0.6, outlier.size = 0.8) +
  labs(x = "Q", y = "sigma2 gap, wb - oracle (% of true)",
       title = "Paired algorithm gap in sigma2")
save_fig(p4c, "step4c_gap_sigma2")
print(p4c)

p4d <- ggplot(pw, aes(Qf, ratio_d)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_boxplot(fill = col_m["wb"], alpha = 0.6, outlier.size = 0.8) +
  labs(x = "Q", y = "d_true(wb) / d_true(oracle)",
       title = "Paired algorithm gap in B")
save_fig(p4d, "step4d_ratio_d")
print(p4d)

step_title("STEP 5. Three-layer decomposition of wb error")

cat("sigma2. total bias = (oracle - truth) + (wb - oracle). Both layers are additive in %.\n")
cat("B. d^2 = 2*sum(sin^2) is additive in first order, so split d_wb^2 = d_oracle^2 + excess.\n\n")

dec_s2 <- pw %>%
  transmute(Qf, seed, `a. truth -> oracle` = sigma2_bias_pct_oracle, `b. oracle -> wb` = gap_s2_pct) %>%
  pivot_longer(-c(Qf, seed), names_to = "layer", values_to = "bias_pct")

dec_B <- pw %>%
  transmute(Qf, seed, `a. truth -> oracle` = d_true_oracle^2, `b. oracle -> wb` = ex_d2) %>%
  pivot_longer(-c(Qf, seed), names_to = "layer", values_to = "d2")

tab_dec <- pw %>%
  group_by(Q) %>%
  summarise(s2_layer_a = mean(sigma2_bias_pct_oracle), s2_layer_b = mean(gap_s2_pct),
            s2_total = mean(sigma2_bias_pct_wb), s2_share_b = s2_layer_b / s2_total,
            d2_layer_a = mean(d_true_oracle^2), d2_layer_b = mean(ex_d2),
            d2_share_b = mean(share_alg), .groups = "drop")
print(tab_dec %>% mutate(across(where(is.double), ~ round(.x, 4))), n = Inf)

p5a <- ggplot(dec_s2, aes(Qf, bias_pct, fill = layer)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c("grey55", "#2c7fb8")) +
  labs(x = "Q", y = "contribution to sigma2 bias (%)", fill = NULL,
       title = "sigma2 error by layer",
       subtitle = "layer a is flat in Q, layer b carries all the Q dependence")
save_fig(p5a, "step5a_layers_sigma2")
print(p5a)

p5b <- ggplot(dec_B, aes(Qf, d2, fill = layer)) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c("grey55", "#2c7fb8")) +
  labs(x = "Q", y = "contribution to d^2", fill = NULL, title = "B error by layer")
save_fig(p5b, "step5b_layers_B")
print(p5b)

p5c <- ggplot(pw, aes(Qf, share_alg)) +
  geom_boxplot(fill = col_m["wb"], alpha = 0.6, outlier.size = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Q", y = "share of wb d^2 from the algorithm layer",
       title = "How much of the wb B error is algorithmic")
save_fig(p5c, "step5c_share_alg")
print(p5c)

step_title("STEP 6. Shape of the Q trend (single J, descriptive only)")

forms <- list(linear_Q = ~ Q, log_Q = ~ log(Q), sqrt_Q = ~ sqrt(Q))
fit_forms <- function(resp) {
  imap_dfr(forms, function(f, nm) {
    m <- lm(update(f, paste(resp, "~ .")), data = pw)
    tibble(response = resp, form = nm, slope = coef(m)[2], AIC = AIC(m),
           resid_cell_means = max(abs(tapply(resid(m), pw$Q, mean))))
  }) %>% mutate(dAIC = AIC - min(AIC))
}
tab_form <- bind_rows(fit_forms("gap_s2_pct"), fit_forms("ratio_d"))
print(tab_form %>% mutate(across(where(is.double), ~ signif(.x, 4))), n = Inf)
cat("Six Q values at one J cannot separate Q from Q/J. Treat the best form as a working description.\n")

best_s2 <- lm(gap_s2_pct ~ log(Q), data = pw)
best_d <- lm(ratio_d ~ Q, data = pw)
fit_pts <- tibble(Q = Q_vals) %>%
  mutate(Qf = factor(Q, levels = Q_vals),
         gap_fit = predict(best_s2, newdata = .),
         ratio_fit = predict(best_d, newdata = .))

p6a <- ggplot(pw, aes(Qf, gap_s2_pct)) +
  geom_boxplot(fill = col_m["wb"], alpha = 0.4, outlier.size = 0.8) +
  geom_line(data = fit_pts, aes(Qf, gap_fit, group = 1), colour = "firebrick") +
  geom_point(data = fit_pts, aes(Qf, gap_fit), colour = "firebrick", size = 2) +
  labs(x = "Q", y = "sigma2 gap (%)", title = "sigma2 gap with a + b log(Q) fit")
save_fig(p6a, "step6a_gap_trend")
print(p6a)

p6b <- ggplot(pw, aes(Qf, ratio_d)) +
  geom_boxplot(fill = col_m["wb"], alpha = 0.4, outlier.size = 0.8) +
  geom_line(data = fit_pts, aes(Qf, ratio_fit, group = 1), colour = "firebrick") +
  geom_point(data = fit_pts, aes(Qf, ratio_fit), colour = "firebrick", size = 2) +
  labs(x = "Q", y = "d ratio wb / oracle", title = "B ratio with a + b Q fit")
save_fig(p6b, "step6b_ratio_trend")
print(p6b)

step_title("STEP 7. Does the algorithm gap depend on signal strength?")

pw <- pw %>%
  group_by(Q) %>%
  mutate(omega_bin = cut(omega_min, quantile(omega_min, c(0, 1 / 3, 2 / 3, 1)),
                         include.lowest = TRUE, labels = c("low", "mid", "high"))) %>%
  ungroup()

tab_om <- pw %>%
  group_by(Q) %>%
  summarise(cor_ratio_logomega = cor(ratio_d, log(omega_min)),
            cor_gap_logomega = cor(gap_s2_pct, log(omega_min)),
            cor_ratio_Ntotal = cor(ratio_d, N_total),
            cor_gap_Ntotal = cor(gap_s2_pct, N_total), .groups = "drop")
print(tab_om %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf)
cat("N_total is the total number of observations in the seed. It measures per-group information.\n")

p7a <- ggplot(pw, aes(Qf, ratio_d, fill = omega_bin)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_boxplot(outlier.size = 0.6, position = position_dodge(0.8), width = 0.7) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Q", y = "d ratio wb / oracle", fill = "omega_min tercile",
       title = "B gap by signal strength within each Q")
save_fig(p7a, "step7a_ratio_by_omega")
print(p7a)

p7b <- ggplot(pw, aes(Qf, gap_s2_pct, fill = omega_bin)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.size = 0.6, position = position_dodge(0.8), width = 0.7) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Q", y = "sigma2 gap (%)", fill = "omega_min tercile",
       title = "sigma2 gap by signal strength within each Q")
save_fig(p7b, "step7b_gap_by_omega")
print(p7b)

step_title("STEP 8. Seed-level outliers (robust z within Q)")

outl <- pw %>%
  group_by(Q) %>%
  mutate(z_ratio = rob_z(ratio_d), z_gap = rob_z(gap_s2_pct)) %>%
  ungroup() %>%
  filter(abs(z_ratio) > 4 | abs(z_gap) > 4) %>%
  transmute(tag = make_tag(J0, Q, "wb"), seed, ratio_d, gap_s2_pct, omega_min, N_total,
            z_ratio = round(z_ratio, 2), z_gap = round(z_gap, 2))
cat(sprintf("%d seeds with |robust z| > 4.\n", nrow(outl)))
print(outl, n = Inf)


step_title("STEP 9. Secondary metrics")

sec <- rr2 %>%
  select(Qf, method, ang_max, tucker, rv, rmse, sig_err) %>%
  pivot_longer(-c(Qf, method), names_to = "metric", values_to = "value")
print(sec %>% group_by(metric, method) %>%
        summarise(min = min(value), median = median(value), max = max(value), .groups = "drop") %>%
        mutate(across(where(is.double), ~ signif(.x, 4))), n = Inf)
cat("tucker and rv sit near 1 for both methods. They separate the methods weakly.\n")
cat("d_true and ang_max are the primary B metrics.\n")

p9 <- ggplot(sec, aes(Qf, value, fill = method)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8), width = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = col_m) +
  labs(x = "Q", y = NULL, fill = NULL, title = "Secondary metrics")
save_fig(p9, "step9_secondary_metrics", w = 11, h = 6.5)
print(p9)

step_title("STEP 10. Convergence and cost caveat for wb")

wbc <- rr2 %>%
  filter(method == "wb") %>%
  mutate(last_step = sigma2_last - sigma2) %>%
  left_join(pw %>% select(Q, seed, gap_s2_pct), by = c("Q", "seed"))

cat(sprintf("wb converged in %.1f%% of runs.\n", 100 * mean(wbc$converged)))
print(wbc %>% group_by(Q) %>%
        summarise(n_max_iter = sum(stop_reason == "max_iter", na.rm = TRUE),
                  n_param_stable = sum(stop_reason == "param_stable", na.rm = TRUE),
                  mean_iter = mean(iterations),
                  frac_peak_is_last = mean(best_iter == iterations), .groups = "drop"), n = Inf)
cat("frac_peak_is_last = 1 means log-evidence was still rising when the run stopped.\n")
cat("The peak snapshot therefore equals the final iterate and cannot help here.\n\n")
cat("sigma2 and sigma2_last differ by one recorded update. Its sign depends on where the\n")
cat("snapshot is taken, so check it against sigma2_trace before reading the direction.\n")
cat("If sigma2 contracts geometrically with ratio r, the remaining drift is about |step|*r/(1-r).\n\n")

left <- wbc %>%
  select(Qf, seed, last_step, gap_s2_pct) %>%
  crossing(r = c(0.97, 0.98, 0.99)) %>%
  mutate(leftover_pct = 100 * abs(last_step) * r / (1 - r) / s2_true,
         leftover_over_gap = leftover_pct / abs(gap_s2_pct))
print(left %>% group_by(Qf, r) %>%
        summarise(med_leftover_pct = median(leftover_pct),
                  med_leftover_over_gap = median(leftover_over_gap), .groups = "drop") %>%
        pivot_wider(names_from = r, values_from = c(med_leftover_pct, med_leftover_over_gap)) %>%
        mutate(across(where(is.double), ~ signif(.x, 3))), n = Inf)
cat("A ratio near or above 1 means the Q=small gap cannot be read before Aitken extrapolation.\n")

p10a <- ggplot(left, aes(Qf, leftover_over_gap, fill = factor(r))) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8), width = 0.7) +
  scale_y_log10() +
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Q", y = "possible leftover drift / |sigma2 gap|", fill = "assumed r",
       title = "Could non-convergence explain the sigma2 gap?")
save_fig(p10a, "step10a_leftover_vs_gap")
print(p10a)

tm <- lm(log(elapsed_sec) ~ log(Q), data = wbc)
cat(sprintf("wb runtime scales as Q^%.2f (100 iterations, J = %d).\n", coef(tm)[2], as.integer(J0)))
p10b <- ggplot(wbc, aes(Qf, elapsed_sec)) +
  geom_boxplot(fill = col_m["wb"], alpha = 0.6, outlier.size = 0.8) +
  scale_y_log10() +
  labs(x = "Q", y = "wall time (s, log scale)",
       title = sprintf("wb runtime, slope %.2f on log-log", coef(tm)[2]))
save_fig(p10b, "step10b_runtime")
print(p10b)

step_title("STEP 11. Summary table")

summary_tab <- tab_or_s2 %>%
  select(Q, oracle_s2_bias = obs_bias, oracle_s2_theory = pred_bias) %>%
  left_join(tab_sig %>% select(Q, omega_true, real_over_true = ratio_mean), by = "Q") %>%
  left_join(tab_or_B %>% select(Q, oracle_d = d_obs, oracle_d_theory = d_pred,
                                oracle_d_theory_true = d_pred_true), by = "Q") %>%
  left_join(tab_pair %>% select(Q, gap_s2_pct, ratio_d, frac_wb_d_worse), by = "Q") %>%
  left_join(tab_dec %>% select(Q, s2_share_b, d2_share_b), by = "Q") %>%
  mutate(across(where(is.double), ~ round(.x, 4)))
print(summary_tab, n = Inf)
write.csv(summary_tab, file.path(fig_dir, sprintf("J%d_summary.csv", as.integer(J0))), row.names = FALSE)
cat("Figures and summary saved in", fig_dir, "\n")
