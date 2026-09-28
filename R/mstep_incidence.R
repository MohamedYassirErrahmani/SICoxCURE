# =============================================================
# mstep_incidence.R
# M-step for the incidence part
# =============================================================

extract_free_params <- function(gamma, constraint) {
  if (constraint == "gamma1") {
    return(gamma[-1L])
  } else {
    return(gamma[-1L])
  }
}

reconstruct_gamma <- function(free, constraint, sign_gamma1 = 1L) {
  if (constraint == "gamma1") {
    gamma <- c(1.0, free)
  } else {
    gamma <- c(sign_gamma1, free)
    nrm   <- sqrt(sum(gamma^2))
    if (nrm > 1e-12) gamma <- gamma / nrm
  }
  gamma
}

neg_loglik_incidence <- function(free, X, W, h, kernel_fn,
                                 constraint, sign_gamma1 = 1L) {
  gamma  <- reconstruct_gamma(free, constraint, sign_gamma1)
  u      <- compute_index(X, gamma)
  g_hat  <- nw_loo(u, W, h, kernel_fn)
  g_hat  <- clip_prob(g_hat)
  -sum(W * safe_log(g_hat) + (1.0 - W) * safe_log(1.0 - g_hat))
}

mstep_incidence <- function(gamma_cur, X, W, h, kernel_fn,
                             constraint) {
  sign_g1 <- if (gamma_cur[1L] >= 0) 1L else -1L
  free0   <- extract_free_params(gamma_cur, constraint)

  # MODIFICATION 1 : supprimer les warnings L-BFGS-B
  # -> gérer en interne avec tryCatch
  opt <- tryCatch({
    stats::optim(
      par     = free0,
      fn      = neg_loglik_incidence,
      method  = "L-BFGS-B",
      X       = X, W = W, h = h,
      kernel_fn = kernel_fn,
      constraint = constraint,
      sign_gamma1 = sign_g1,
      control = list(maxit = 200L, factr = 1e10)
    )
  }, warning = function(w) {
    # L-BFGS-B warning -> try Nelder-Mead silently
    stats::optim(
      par    = free0,
      fn     = neg_loglik_incidence,
      method = "Nelder-Mead",
      X = X, W = W, h = h,
      kernel_fn = kernel_fn,
      constraint = constraint,
      sign_gamma1 = sign_g1,
      control = list(maxit = 500L)
    )
  }, error = function(e) {
    # If both fail, return current parameters
    list(par = free0, convergence = 1L)
  })

  # If L-BFGS-B didn't converge (code != 0), try Nelder-Mead
  if (opt$convergence != 0L) {
    opt <- tryCatch({
      stats::optim(
        par    = free0,
        fn     = neg_loglik_incidence,
        method = "Nelder-Mead",
        X = X, W = W, h = h,
        kernel_fn = kernel_fn,
        constraint = constraint,
        sign_gamma1 = sign_g1,
        control = list(maxit = 500L)
      )
    }, error = function(e) {
      list(par = free0, convergence = 1L)
    })
  }

  gamma_new <- reconstruct_gamma(opt$par, constraint, sign_g1)
  u_new     <- compute_index(X, gamma_new)
  g_hat     <- nw_loo(u_new, W, h, kernel_fn)
  g_hat     <- clip_prob(g_hat)

  list(gamma = gamma_new, g_hat = g_hat)
}
