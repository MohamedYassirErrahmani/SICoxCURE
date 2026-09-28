# =============================================================
# mstep_latency.R
# M-step for the latency : weighted Cox + Breslow
# =============================================================

# ── Score to MINIMISE for a given beta ───────────────────────
neg_partial_loglik <- function(beta, Y, delta, Z, W) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)

  lp         <- as.vector(Z %*% beta)

  # Numerical stabilisation : subtract max(lp) before exp()
  # Does not change the score (exp(-max) cancels in num/denom)
  lp         <- lp - max(lp)

  # Weighted risk score : W_k * exp(beta^t * Z_k)
  risk_score <- W *safe_exp(lp)

  # Sort by Y decreasing to compute cumulative risk via cumsum()
  # cum_risk[i] = sum_{k in Ri} W_k * exp(beta^t * Z_k)
  n              <- length(Y)
  ord            <- order(Y, decreasing = TRUE)
  cum_risk       <- numeric(n)
  cum_risk[ord]  <- cumsum(risk_score[ord])
  log_cum_risk   <-safe_log(cum_risk)

  # log(L2_breve) = sum_i delta_i * (lp_i - log(cum_risk_i))
  # delta acts as switch : delta=1 -> contributes, delta=0 -> 0
  ll <- sum(delta * (lp - log_cum_risk))

  # Return -log(L2_breve) : minimising this = maximising L2_breve
  -ll
}

# ── Main function — find the best beta ───────────────────────
mstep_latency <- function(beta_cur, Y, delta, Z, W,
                          max_iter = 200L) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)

  # ── Primary optimisation : L-BFGS-B ─────────────────────────
  # numerical techniques such as the Newton-Raphson algorithm"
  opt <- tryCatch(
    stats::optim(
      par     = beta_cur,
      fn      = neg_partial_loglik,
      Y       = Y, delta = delta, Z = Z, W = W,
      method  = "L-BFGS-B",
      control = list(maxit = max_iter)
    ),
    error = function(e) NULL
  )

  # ── Fallback : Nelder-Mead if L-BFGS-B failed ───────────────
  if (is.null(opt) || opt$convergence != 0L) {

    if (is.null(opt)) {
      warning(sprintf(
        "mstep_latency: L-BFGS-B optimisation failed completely.
         Trying Nelder-Mead as fallback (max_iter=%d).",
        max_iter * 2L
      ))
    } else {
      warning(sprintf(
        "mstep_latency: L-BFGS-B did not converge
         (convergence code=%d). Trying Nelder-Mead as fallback.",
        opt$convergence
      ))
    }

    opt <- tryCatch(
      stats::optim(
        par     = beta_cur,
        fn      = neg_partial_loglik,
        Y       = Y, delta = delta, Z = Z, W = W,
        method  = "Nelder-Mead",
        control = list(maxit = max_iter * 2L)
      ),
      error = function(e) NULL
    )
  }

  # ── If both methods failed : keep current beta ───────────────
  if (!is.null(opt)) {
    beta_new <- opt$par
    ll_lat   <- -opt$value
  } else {
    warning(
      "mstep_latency: both L-BFGS-B and Nelder-Mead failed.
       Keeping current beta unchanged."
    )
    beta_new <- beta_cur
    ll_lat   <- NA_real_
  }

  # ── Recalculate breslow and Su_hat with the best beta ────────
  breslow <- breslow_estimator(Y, delta, Z, beta_new, W)

  # Su_hat_i = Su(Yi|Zi) for each patient, used by estep()
  # at the next iteration to recalculate W
  Su_hat  <- pmax(eval_Su_diag(Y, Z, beta_new, breslow), 0.0)

  list(beta    = beta_new,
       breslow = breslow,
       Su_hat  = Su_hat,
       ll_lat  = ll_lat)
}
