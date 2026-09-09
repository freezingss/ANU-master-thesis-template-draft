source("basic_functions.R")

fhem_loglik <- function(Sigma, lt, Fl_inv) {
  total <- 0
  for (j in seq_along(lt)) {
    V <- Sigma + Fl_inv[[j]]
    ldet <- as.numeric(determinant(V, logarithm = TRUE)$modulus)
    total <- total - 0.5 * ldet - 0.5 * drop(crossprod(lt[[j]], solve(V, lt[[j]])))
  }
  total
}

fhem_backsolve <- function(lambda_hat, S_hat, tau2, ridge = 1e-3) {
  Q <- nrow(lambda_hat)
  J <- ncol(lambda_hat)
  P <- diag(Q) / tau2
  lt <- vector("list", J)
  Fl <- vector("list", J)
  Fl_inv <- vector("list", J)
  nbad <- 0L
  for (j in 1:J) {
    Sh <- S_hat[[j]]
    res <- tryCatch({
      Sh_inv <- solve(Sh)
      Fj <- Sh_inv - P
      Fj_reg <- Fj + diag(ridge, Q)
      list(lt = solve(Fj_reg, Sh_inv %*% lambda_hat[, j]), Fl = Fj_reg, Fl_inv = solve(Fj_reg))
    }, error = function(e) NULL)
    if (is.null(res)) {
      nbad <- nbad + 1L
      Fl[[j]] <- diag(ridge, Q)
      Fl_inv[[j]] <- diag(Q) / ridge
      lt[[j]] <- rep(0, Q)
    } else {
      lt[[j]] <- res$lt
      Fl[[j]] <- res$Fl
      Fl_inv[[j]] <- res$Fl_inv
    }
  }
  if (nbad > 0) warning(sprintf("fhem_backsolve: %d of %d groups failed (singular S_hat) and were treated as uninformative", nbad, J))
  list(lt = lt, Fl = Fl, Fl_inv = Fl_inv, n_failed = nbad)
}

fhem_gaussian_em <- function(lt, Fl, Fl_inv, Sigma0, K, max_iter = 300, tol = 1e-6) {
  J <- length(lt)
  q <- nrow(Sigma0)
  Sigma <- Sigma0
  ll <- fhem_loglik(Sigma, lt, Fl_inv)
  trace <- ll
  fit <- ppca_closed(Sigma, K)
  it <- 0
  for (it in 1:max_iter) {
    Sinv <- tryCatch(solve(Sigma), error = function(e) {
      warning(sprintf("fhem_gaussian_em: solve(Sigma) failed at iter %d (%s); falling back to a diagonal approximation", it, conditionMessage(e)))
      diag(1 / pmax(diag(Sigma), 1e-8))
    })
    S <- matrix(0, q, q)
    for (j in 1:J) {
      V <- solve(Fl[[j]] + Sinv)
      m <- V %*% (Fl[[j]] %*% lt[[j]])
      S <- S + tcrossprod(m) + V
    }
    S <- S / J
    fit <- ppca_closed(S, K)
    Sigma <- fit$Sigma
    ll_new <- fhem_loglik(Sigma, lt, Fl_inv)
    trace <- c(trace, ll_new)
    if (abs(ll_new - ll) < tol) break
    ll <- ll_new
  }
  list(fit = fit, trace = trace, iters = it, monotone = all(diff(trace) > -1e-6))
}

fhem_from_iso_estep <- function(lambda_hat, S_hat, tau2, K,
                                ridge = 1e-3, max_iter = 300, tol = 1e-6) {
  Q <- nrow(lambda_hat)
  bs <- fhem_backsolve(lambda_hat, S_hat, tau2, ridge = ridge)
  em <- fhem_gaussian_em(bs$lt, bs$Fl, bs$Fl_inv, diag(Q) * tau2, K, max_iter = max_iter, tol = tol)
  list(B = em$fit$B, sigma2 = em$fit$sigma2, Sigma = em$fit$Sigma,
       lt = bs$lt, Fl = bs$Fl, trace = em$trace, iters = em$iters, monotone = em$monotone,
       n_failed_groups = bs$n_failed)
}

fit_fhem <- function(K, estep_iso_fn, mstep_phi_fn, apply_PLT_fn = identity,
                     mu0, phi0, B0, sigma2_0,
                     ridge = 1e-3, max_iter = 30, tol_outer = 1e-6,
                     fhem_max_iter = 300, fhem_tol = 1e-6, verbose = TRUE,
                     trace = TRUE) {
  Q <- length(mu0)
  mu <- mu0
  phi <- phi0
  B <- B0
  sigma2 <- sigma2_0
  tau2 <- sum(diag(tcrossprod(B))) / Q + sigma2
  hist_df <- data.frame(iter = integer(0), tau2 = numeric(0), sigma2 = numeric(0), fh_ll = numeric(0))
  fh <- NULL
  
  if (trace) {
    best_iter <- NA_integer_
    best_ll <- -Inf
    best_state <- NULL
  }
  
  for (outer in 1:max_iter) {
    es <- estep_iso_fn(mu, phi, tau2)
    mp <- mstep_phi_fn(es$lambda_hat, mu, phi)
    mu <- mp$mu
    phi <- mp$phi
    fh <- fhem_from_iso_estep(es$lambda_hat, es$S_hat, tau2, K,
                              ridge = ridge, max_iter = fhem_max_iter, tol = fhem_tol)
    B <- apply_PLT_fn(fh$B)
    sigma2 <- fh$sigma2
    tau2_new <- sum(diag(tcrossprod(B))) / Q + sigma2
    ll <- tail(fh$trace, 1)
    hist_df <- rbind(hist_df, data.frame(iter = outer, tau2 = tau2, sigma2 = sigma2, fh_ll = ll))
    
    if (trace && ll > best_ll) {
      best_ll <- ll
      best_iter <- outer
      best_state <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2, tau2 = tau2_new,
                         lambda_hat = es$lambda_hat, S_hat = es$S_hat, lt = fh$lt, Fl = fh$Fl)
    }
    
    if (verbose) cat(sprintf(
      "outer %2d  tau2=%.4f  sigma2=%.4f  fh_ll=%.2f  inner_iters=%d  inner_monotone=%s\n",
      outer, tau2, sigma2, ll, fh$iters, fh$monotone))
    converged <- outer > 1 && abs(tau2_new - tau2) < tol_outer
    tau2 <- tau2_new
    if (converged) break
  }
  
  es_final <- estep_iso_fn(mu, phi, tau2)
  fh_final <- fhem_from_iso_estep(es_final$lambda_hat, es_final$S_hat, tau2, K,
                                  ridge = ridge, max_iter = fhem_max_iter, tol = fhem_tol)
  B <- apply_PLT_fn(fh_final$B)
  sigma2 <- fh_final$sigma2
  tau2 <- sum(diag(tcrossprod(B))) / Q + sigma2
  fh <- fh_final
  fh_ll_final <- tail(fh_final$trace, 1)
  
  out <- list(mu = mu, phi = phi, B = B, sigma2 = sigma2, tau2 = tau2,
              lt = fh$lt, Fl = fh$Fl, trace = hist_df, fh_ll_final = fh_ll_final)
  
  if (trace) {
    out <- c(out, list(best_iter = best_iter, best_ll = best_ll, best = best_state))
  }
  out
}