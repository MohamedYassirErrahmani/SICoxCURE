# =============================================================
# em_algorithm.R
# Main EM algorithm loop
# =============================================================

run_em <- function(Y, delta, X, Z,
                   init, data_info,
                   kernel      = "epanechnikov",
                   h_range     = c(0.4, 1.0),
                   h_step      = 0.001,
                   constraint  = "gamma1",
                   tol         = 1e-5,
                   max_iter    = 500L,
                   conv_method = c("both", "max", "sum"),
                   verbose     = 1L) {

  # ── Initialisation ────────────────────────────────────────────
  conv_method <- match.arg(conv_method)
  kernel_fn   <- get_kernel(kernel)
  beyond_tau  <- data_info$beyond_tau
  n           <- length(Y)

  gamma_cur <- init$gamma
  beta_cur  <- init$beta
  p_cur     <- init$p_hat
  Su_cur    <- init$Su_hat
  h_cur     <- init$h
  breslow   <- init$breslow

  # ── Diagnostic vector ─────────────────────────────────────────
  loglik_comp_vec <- rep(NA_real_, max_iter)

  converged <- FALSE
  n_iter    <- max_iter
  crit      <- Inf
  S0_prev   <- breslow$S0

  if (verbose >= 1L) {
    cat("--- Starting EM algorithm ---\n")
    cat(sprintf("Convergence criterion : %s\n", conv_method))
    cat(sprintf("Bandwidth grid        : [%.2f, %.2f] step=%.4f (%d values)\n",
                h_range[1], h_range[2], h_step,
                length(seq(h_range[1], h_range[2], by=h_step))))

    # MODIFICATION 2 : barre de progression toujours à 100%
    # On utilise max = max_iter mais on la complète à la fin
    pb <- utils::txtProgressBar(min=0, max=100L, style=3)
  }

  for (m in seq_len(max_iter)) {

    # ── E STEP ───────────────────────────────────────────────────
    W <- estep(p_cur, Su_cur, delta, beyond_tau)

    # ── M STEP 1 : Bandwidth ─────────────────────────────────────
    u_cur  <- compute_index(X, gamma_cur)
    cv_res <- select_bandwidth(u_cur, W, kernel_fn,
                               h_range, h_step)
    h_new  <- cv_res$h_opt

    # ── M STEP 2 : Incidence ─────────────────────────────────────
    # Suppress L-BFGS-B warning -> handled internally
    inc_res   <- suppressWarnings(
      mstep_incidence(gamma_cur, X, W, h_new,
                      kernel_fn, constraint)
    )
    gamma_new <- inc_res$gamma
    g_new     <- inc_res$g_hat

    # ── M STEP 3 : Latency ───────────────────────────────────────
    lat_res  <- mstep_latency(beta_cur, Y, delta, Z, W)
    beta_new <- lat_res$beta
    Su_new   <- lat_res$Su_hat
    breslow  <- lat_res$breslow

    # ── CONVERGENCE CRITERION ─────────────────────────────────────
    S0_new      <- breslow$S0
    conv_result <- convergence_criterion(
      beta_new, beta_cur,
      gamma_new, gamma_cur,
      S0_new, S0_prev,
      method = conv_method
    )
    crit <- conv_result$crit

    # ── DIAGNOSTIC ───────────────────────────────────────────────
    hu_approx          <- compute_hazard_u(Y, delta, Z,
                                           beta_new, breslow)
    loglik_comp_vec[m] <- complete_loglik(g_new, Su_new,
                                          hu_approx, delta, W)

    # ── Verbose detail if verbose >= 2 ───────────────────────────
    if (verbose >= 2L) {
      ll_comp_str <- if (!is.na(loglik_comp_vec[m]))
        sprintf("L_c=%.3f", loglik_comp_vec[m])
      else "L_c=NA"
      cat(sprintf(
        "\n  iter=%d | crit=%.2e (%s) | %s | h=%.4f",
        m, crit, conv_result$method_used, ll_comp_str, h_new
      ))
    }

    # ── Update parameters ─────────────────────────────────────────
    S0_prev   <- S0_new
    gamma_cur <- gamma_new
    beta_cur  <- beta_new
    p_cur     <- g_new
    Su_cur    <- Su_new
    h_cur     <- h_new

    # MODIFICATION 2 : progression proportionnelle à la convergence
    # -> représente à quel point on est proche de tol
    if (verbose >= 1L) {
      # Progression = log(crit_init/crit) / log(crit_init/tol)
      # Toujours entre 0 et 100
      if (m == 1L) crit_init <- max(crit, tol * 10)
      pct <- min(99L, as.integer(
        100 * log(crit_init / max(crit, tol * 0.01)) /
              log(crit_init / tol)
      ))
      pct <- max(0L, pct)
      utils::setTxtProgressBar(pb, pct)
    }

    # ── Convergence test ──────────────────────────────────────────
    if (crit < tol && m > 1L) {
      converged <- TRUE
      n_iter    <- m
      break
    }
  }

  # ── End of loop ───────────────────────────────────────────────
  if (verbose >= 1L) {
    # Barre à 100% quand convergé
    utils::setTxtProgressBar(pb, 100L)
    close(pb)
    cat("\n")

    if (converged) {
      # MODIFICATION 1 : message en vert + clair
      green <- "\033[32m"
      reset <- "\033[0m"
      cat(sprintf(
        "%s\u2713 Model converged successfully at iteration %d (criterion = %.2e < tol = %.0e)%s\n",
        green, n_iter, crit, tol, reset
      ))
    } else {
      warning(sprintf(
        "No convergence after %d iterations (crit=%.2e).",
        max_iter, crit
      ))
    }

    if (!is.na(loglik_comp_vec[n_iter]))
      cat(sprintf("  Final complete log-likelihood : %.4f\n",
                  loglik_comp_vec[n_iter]))
  }

  # ── Return results ────────────────────────────────────────────
  list(
    gamma       = gamma_cur,
    beta        = beta_cur,
    g_hat       = p_cur,
    p_hat       = p_cur,
    Su_hat      = Su_cur,
    breslow     = breslow,
    h           = h_cur,
    W           = W,
    loglik      = loglik_comp_vec[seq_len(n_iter)],
    loglik_comp = loglik_comp_vec[seq_len(n_iter)],
    converged   = converged,
    n_iter      = n_iter,
    crit        = crit,
    conv_method = conv_method
  )
}
