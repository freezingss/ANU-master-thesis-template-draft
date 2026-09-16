source("basic_functions.R")
source("pfa_debias.R")

set.seed(42)

Q <- 5
P <- 3
K <- 2
J <- 300
Nj <- 30
N <- J * Nj

mu_true <- c(0, 0.3, -0.2, 0.5, -0.1)
phi_true <- rbind(
  rep(0, Q),
  c(0, 0.8, -0.6, 1.2, -1.0),
  c(0, -0.5, 0.9, -1.3, 0.7)
)
B_true <- matrix(c(0.3, -0.2, 0.1, 0.4, -0.1,
                   0.2, 0.3, -0.4, 0.1, 0.2), nrow = Q, ncol = K)
sigma2_true <- 0.15

group <- rep(1:J, each = Nj)
X <- cbind(1, rnorm(N, 0, 1), rnorm(N, 0, 1))

lambda <- matrix(0, Q, J)
for (j in 1:J) {
  f_j <- rnorm(K)
  eps_j <- rnorm(Q, 0, sqrt(sigma2_true))
  lambda[, j] <- B_true %*% f_j + eps_j
}

M <- sample(200:300, N, replace = TRUE)
Y <- matrix(0, N, Q)
for (i in 1:N) {
  eta_i <- mu_true + as.numeric(X[i, ] %*% phi_true) + lambda[, group[i]]
  eta_i <- eta_i - max(eta_i)
  p_i <- exp(eta_i) / sum(exp(eta_i))
  Y[i, ] <- as.numeric(rmultinom(1, M[i], p_i))
}

long <- build_long(Y, X, group)
xnames <- paste0("x", 2:P)
covar_terms <- paste(paste0("category:", xnames), collapse = " + ")
form <- as.formula(paste0("count ~ category + ", covar_terms,
                          " + rr(category + 0 | group, d = ", K, ") + (1 | obs)"))

fit <- glmmTMB(form, data = long, family = poisson(link = "log"),
               offset = long$log_total,
               control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

cat("converged (pdHess):", isTRUE(fit$sdr$pdHess), "\n\n")

fe <- fixef(fit)$cond
print(fe)

extract_by_cat <- function(fe, xname, Q) {
  pat <- paste0("^category([0-9]+):", xname, "$")
  nm <- names(fe)
  idx <- grep(pat, nm)
  cats <- as.integer(sub(pat, "\\1", nm[idx]))
  out <- rep(NA_real_, Q)
  out[cats] <- as.numeric(fe[idx])
  out
}

phi2_hat <- extract_by_cat(fe, "x2", Q)
phi3_hat <- extract_by_cat(fe, "x3", Q)

cat("\nx2 slopes, true vs estimated:\n")
print(data.frame(category = 1:Q, true = phi_true[2, ], estimated = phi2_hat))
cat("\nx3 slopes, true vs estimated:\n")
print(data.frame(category = 1:Q, true = phi_true[3, ], estimated = phi3_hat))

cat("\nnumber of x2 slopes recovered:", sum(!is.na(phi2_hat)), " (expected Q =", Q, ")\n")
cat("number of x3 slopes recovered:", sum(!is.na(phi3_hat)), " (expected Q =", Q, ")\n")
cat("correlation x2 (true vs estimated):", cor(phi_true[2, ], phi2_hat, use = "complete.obs"), "\n")
cat("correlation x3 (true vs estimated):", cor(phi_true[3, ], phi3_hat, use = "complete.obs"), "\n")

# Remove rr()
form_no_rr <- as.formula(paste0("count ~ category + ", covar_terms, " + (1 | obs)"))

fit_no_rr <- glmmTMB(form_no_rr, data = long, family = poisson(link = "log"),
                     offset = long$log_total,
                     control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

fe2 <- fixef(fit_no_rr)$cond
phi2_hat_norr <- extract_by_cat(fe2, "x2", Q)
phi3_hat_norr <- extract_by_cat(fe2, "x3", Q)

fit_line <- function(t, e) {
  A <- cbind(t, 1)
  coef(lm.fit(A, e))
}
cat("x2 slope without rr():", fit_line(phi_true[2, ], phi2_hat_norr)[1], "\n")
cat("x3 slope without rr():", fit_line(phi_true[3, ], phi3_hat_norr)[1], "\n")

form_fixed_only <- as.formula(paste0("count ~ category + ", covar_terms))
fit_fixed_only <- glm(form_fixed_only, data = long, family = poisson(link = "log"),
                      offset = long$log_total)

fe3 <- coef(fit_fixed_only)
phi2_hat_fixedonly <- extract_by_cat(fe3, "x2", Q)
phi3_hat_fixedonly <- extract_by_cat(fe3, "x3", Q)

cat("x2 slope, fixed-effects-only glm:", fit_line(phi_true[2, ], phi2_hat_fixedonly)[1], "\n")
cat("x3 slope, fixed-effects-only glm:", fit_line(phi_true[3, ], phi3_hat_fixedonly)[1], "\n")

library(nnet)

df_wide <- data.frame(x2 = X[, 2], x3 = X[, 3])
fit_multinom <- multinom(Y ~ x2 + x3, data = data.frame(x2 = X[, 2], x3 = X[, 3]), trace = FALSE)

co <- coef(fit_multinom)
print(co)
phi2_true_contrast <- phi_true[2, ] - phi_true[2, 1]
phi3_true_contrast <- phi_true[3, ] - phi_true[3, 1]

phi2_hat_multinom <- c(0, co[, "x2"])
phi3_hat_multinom <- c(0, co[, "x3"])

cat("x2 slope, nnet::multinom:", fit_line(phi2_true_contrast, phi2_hat_multinom)[1], "\n")
cat("x3 slope, nnet::multinom:", fit_line(phi3_true_contrast, phi3_hat_multinom)[1], "\n")