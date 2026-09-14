# =============================================================
# mstep_incidence.R
# M-step for the incidence : optimisation of gamma
# =============================================================

# ── STEP 1 : Extract free parameters from gamma ──────────────
# Before optim() : remove the fixed parameter (gamma[1] if
# constraint="gamma1") so optim() only optimises the free ones
extract_free_params <- function(gamma, constraint) {
  if (constraint == "gamma1") gamma[-1L] else gamma
}

# ── STEP 2 : Reconstruct full gamma from free parameters ─────
# After optim() : put gamma[1] back in its place
reconstruct_gamma <- function(gamma_free, constraint, d) {
  if (constraint == "gamma1") {
    c(1.0, gamma_free)
  } else {
    nrm <- sqrt(sum(gamma_free^2))
    if (nrm < 1e-12) {
      warning(sprintf(
        "reconstruct_gamma: norm of gamma_free is nearly zero (%.2e).
         Returning default gamma = c(1, 0, ..., 0).",
        nrm
      ))
      return(c(1.0, rep(0.0, d - 1L)))
    }
    g <- gamma_free / nrm
    if (g[1] < 0.0) g <- -g
    g
  }
}

# ── STEP 3 : Score to MINIMISE for a given gamma_free ────────
neg_loglik_incidence <- function(gamma_free, X, W, h,
                                 kernel_fn, constraint) {
  gamma <- reconstruct_gamma(gamma_free, constraint, ncol(X))
  u     <- compute_index(X, gamma)
  g_hat <- nw_loo(u, W, h, kernel_fn)
  -sum(W *safe_log(g_hat) + (1.0 - W) *safe_log(1.0 - g_hat))
}

# ── STEP 4 : Main function — find the best gamma ─────────────
mstep_incidence <- function(gamma_cur, X, W, h, kernel_fn,
                            constraint   = "gamma1",
                            max_iter_opt = 200L) {
  if (!is.matrix(X)) X <- as.matrix(X)
  d <- ncol(X)

  # Special case d=1 : gamma is fixed to 1 regardless of constraint
  # -> constraint="gamma1" : gamma[1]=1, nothing to optimise
  # -> constraint="norm1"  : ||gamma||=1 and gamma[1]>0 -> gamma=1
  # reduces to a non-parametric estimator"
  if (d == 1L) {
    gamma_new <- 1.0
    u         <- compute_index(X, gamma_new)
    g_hat     <- nw_loo(u, W, h, kernel_fn)
    ll        <- sum(W *safe_log(g_hat) +
                       (1.0 - W) *safe_log(1.0 - g_hat))
    return(list(gamma = gamma_new, g_hat = g_hat, ll_inc = ll))
  }

  # Extract free parameters from current gamma
  gamma_free_init <- extract_free_params(gamma_cur, constraint)

  # ── Primary optimisation : L-BFGS-B ─────────────────────────
  # numerical techniques such as the Newton-Raphson algorithm"
  opt <- tryCatch(
    stats::optim(
      par       = gamma_free_init,
      fn        = neg_loglik_incidence,
      X         = X, W = W, h = h,
      kernel_fn = kernel_fn,
      constraint = constraint,
      method    = "L-BFGS-B",
      control   = list(maxit = max_iter_opt)
    ),
    error = function(e) NULL
  )

  # ── Fallback : Nelder-Mead if L-BFGS-B failed ───────────────
  if (is.null(opt) || opt$convergence != 0L) {

    if (is.null(opt)) {
      warning(sprintf(
        "mstep_incidence: L-BFGS-B optimisation failed completely.
         Trying Nelder-Mead as fallback (max_iter=%d).",
        max_iter_opt * 2L
      ))
    } else {
      warning(sprintf(
        "mstep_incidence: L-BFGS-B did not converge
         (convergence code=%d). Trying Nelder-Mead as fallback.",
        opt$convergence
      ))
    }

    opt <- tryCatch(
      stats::optim(
        par       = gamma_free_init,
        fn        = neg_loglik_incidence,
        X         = X, W = W, h = h,
        kernel_fn = kernel_fn,
        constraint = constraint,
        method    = "Nelder-Mead",
        control   = list(maxit = max_iter_opt * 2L)
      ),
      error = function(e) NULL
    )
  }

  # ── If both methods failed : keep current gamma ──────────────
  if (!is.null(opt)) {
    gamma_new <- reconstruct_gamma(opt$par, constraint, d)
    gamma_new <- normalize_gamma(gamma_new, constraint)$gamma
  } else {
    warning(
      "mstep_incidence: both L-BFGS-B and Nelder-Mead failed.
       Keeping current gamma unchanged."
    )
    gamma_new <- gamma_cur
  }

  # ── Final g_hat with the best gamma found ────────────────────
  # is estimated by the estimator given in (4), but with gamma
  # replaced by gamma_hat^(m)"
  u     <- compute_index(X, gamma_new)
  g_hat <- nw_loo(u, W, h, kernel_fn)
  ll    <- sum(W *safe_log(g_hat) +
                 (1.0 - W) *safe_log(1.0 - g_hat))

  list(gamma = gamma_new, g_hat = g_hat, ll_inc = ll)
}
