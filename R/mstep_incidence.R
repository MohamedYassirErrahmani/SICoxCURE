# =============================================================
# mstep_incidence.R  (_bis version)
# M-step for the incidence part
#
# constraint = "gamma1" :
#   -> fix gamma[1] = 1 or -1
#   -> optimize (gamma2,...,gammad) freely with BFGS
#   -> Nelder-Mead fallback
#
# constraint = "norm1" :
#   -> optimize gamma under ||gamma|| = 1
#   -> nloptr SLSQP (same as final version)
#   -> BFGS with re-normalization fallback
#   -> sign of gamma[1] fixed after optimization
# =============================================================

# ── GAMMA1 : reconstruct gamma ────────────────────────────────
reconstruct_gamma_gamma1 <- function(free, sign_gamma1 = 1L) {
  c(as.numeric(sign_gamma1), free)
}

# ── neg_loglik for gamma1 ─────────────────────────────────────
neg_loglik_gamma1 <- function(free, X, W, h, kernel_fn,
                              sign_gamma1 = 1L) {
  gamma <- reconstruct_gamma_gamma1(free, sign_gamma1)
  u     <- compute_index(X, gamma)
  g_hat <- nw_loo(u, W, h, kernel_fn)
  g_hat <- clip_prob(g_hat)
  -sum(W * safe_log(g_hat) + (1.0 - W) * safe_log(1.0 - g_hat))
}

# ── neg_loglik for norm1 (direct on gamma) ────────────────────
neg_loglik_norm1 <- function(gamma, X, W, h, kernel_fn) {
  u     <- compute_index(X, gamma)
  g_hat <- nw_loo(u, W, h, kernel_fn)
  g_hat <- clip_prob(g_hat)
  -sum(W * safe_log(g_hat) + (1.0 - W) * safe_log(1.0 - g_hat))
}

# ── neg_loglik + gradient for nloptr ─────────────────────────
neg_loglik_norm1_grad <- function(gamma, X, W, h, kernel_fn) {
  list(
    "objective" = neg_loglik_norm1(gamma, X, W, h, kernel_fn),
    "gradient"  = nloptr::nl.grad(
      gamma,
      neg_loglik_norm1,
      X = X, W = W, h = h, kernel_fn = kernel_fn
    )
  )
}

# ── Equality constraint : ||gamma||^2 = 1 ────────────────────
constraint_norm1 <- function(gamma, X, W, h, kernel_fn) {
  list(
    "constraints" = sum(gamma^2) - 1,
    "jacobian"    = 2 * gamma
  )
}

# ── Main M-step ───────────────────────────────────────────────
mstep_incidence <- function(gamma_cur, X, W, h, kernel_fn,
                             constraint) {

  sign_g1 <- if (gamma_cur[1L] >= 0) 1L else -1L

  if (constraint == "gamma1") {

    # ── gamma1 : d-1 free parameters, BFGS ───────────────────
    free0 <- gamma_cur[-1L]

    opt <- tryCatch({
      stats::optim(
        par        = free0,
        fn         = neg_loglik_gamma1,
        method     = "BFGS",
        X          = X, W = W, h = h,
        kernel_fn  = kernel_fn,
        sign_gamma1 = sign_g1,
        control    = list(maxit = 300L, reltol = 1e-10)
      )
    }, warning = function(w) {
      list(par = free0, convergence = 1L)
    }, error = function(e) {
      list(par = free0, convergence = 1L)
    })

    if (is.null(opt) || opt$convergence != 0L) {
      opt <- tryCatch({
        stats::optim(
          par        = free0,
          fn         = neg_loglik_gamma1,
          method     = "Nelder-Mead",
          X          = X, W = W, h = h,
          kernel_fn  = kernel_fn,
          sign_gamma1 = sign_g1,
          control    = list(maxit = 500L)
        )
      }, error = function(e) list(par = free0, convergence = 1L))
    }

    gamma_new <- reconstruct_gamma_gamma1(opt$par, sign_g1)

  } else {

    # ── norm1 : SLSQP sous contrainte ||gamma||=1 ─────────────
    gamma0 <- gamma_cur / sqrt(sum(gamma_cur^2))

    # Tentative 1 : AUGLAG + L-BFGS (sous-optimiseur)
    opt <- tryCatch({
      nloptr::nloptr(
        x0       = gamma0,
        eval_f   = neg_loglik_norm1_grad,
        eval_heq = constraint_norm1,
        opts     = list(
          algorithm   = "NLOPT_LD_AUGLAG",
          xtol_rel    = 1e-8,
          ftol_rel    = 1e-8,
          maxeval     = 1000L,
          print_level = 0L,
          local_opts  = list(
            algorithm = "NLOPT_LD_LBFGS",
            xtol_rel  = 1e-8,
            ftol_rel  = 1e-8,
            maxeval   = 300L
          )
        ),
        X = X, W = W, h = h, kernel_fn = kernel_fn
      )
    }, error = function(e) NULL)

    # Tentative 2 : SLSQP si AUGLAG echoue ou degrade
    if (is.null(opt) || opt$status < 0L) {
      opt <- tryCatch({
        nloptr::nloptr(
          x0       = gamma0,
          eval_f   = neg_loglik_norm1_grad,
          eval_heq = constraint_norm1,
          opts     = list(
            algorithm   = "NLOPT_LD_SLSQP",
            xtol_rel    = 1e-8,
            ftol_rel    = 1e-8,
            maxeval     = 500L,
            print_level = 0L
          ),
          X = X, W = W, h = h, kernel_fn = kernel_fn
        )
      }, error = function(e) NULL)
    }

    if (!is.null(opt) && opt$status >= 0L) {
      gamma_new <- opt$solution
    } else {
      # Fallback : BFGS avec re-normalisation
      opt_fb <- tryCatch({
        stats::optim(
          par     = gamma0,
          fn      = function(g, ...) {
            g <- g / sqrt(sum(g^2))
            neg_loglik_norm1(g, ...)
          },
          method  = "BFGS",
          X = X, W = W, h = h, kernel_fn = kernel_fn,
          control = list(maxit = 300L)
        )
      }, error = function(e) list(par = gamma0))
      gamma_new <- opt_fb$par
    }

    gamma_new <- gamma_new / sqrt(sum(gamma_new^2))
    if (sign(gamma_new[1L]) != sign_g1)
      gamma_new <- -gamma_new
  }

  u_new <- compute_index(X, gamma_new)
  g_hat <- nw_loo(u_new, W, h, kernel_fn)
  g_hat <- clip_prob(g_hat)

  list(gamma = gamma_new, g_hat = g_hat)
}
