# =============================================================
# initialize.R
# EM algorithm initialisation
# =============================================================

init_incidence <- function(X, delta, constraint = "gamma1") {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)
  df        <- as.data.frame(cbind(y = delta, X))
  logit_fit <- tryCatch(
    stats::glm(y ~ ., data = df, family = stats::binomial("logit")),
    error = function(e) NULL
  )
  if (is.null(logit_fit)) {
    gamma_init        <- c(1.0, rep(0.0, d - 1L))
    names(gamma_init) <- colnames(X)
    return(list(gamma = gamma_init, p_hat = rep(mean(delta), n)))
  }
  coefs     <- stats::coef(logit_fit)
  gamma_raw <- coefs[names(coefs) != "(Intercept)"]
  gamma_raw[is.na(gamma_raw)] <- 0.0
  names(gamma_raw) <- colnames(X)
  gamma_init <- normalize_gamma(gamma_raw, constraint)$gamma
  p_hat      <- clip_prob(as.vector(stats::fitted(logit_fit)))
  list(gamma = gamma_init, p_hat = p_hat, logit_fit = logit_fit)
}

init_latency <- function(Y, delta, Z) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  n       <- length(Y)
  z_names <- if (!is.null(colnames(Z))) colnames(Z) else
    paste0("Z", seq_len(ncol(Z)))

  #    observations only for the latency" ──────────────────────
  idx_events <- which(delta == 1L)
  n_events   <- length(idx_events)

  if (n_events < 2L) {
    # Not enough events -> fall back to all observations
    warning(sprintf(
      "Only %d event(s) observed. Using all observations for Cox initialisation.",
      n_events
    ))
    idx_fit <- seq_len(n)
  } else {
    idx_fit <- idx_events
  }

  Y_fit     <- Y[idx_fit]
  delta_fit <- delta[idx_fit]
  Z_fit     <- Z[idx_fit, , drop=FALSE]

  # Build data.frame for coxph()
  df_cox <- as.data.frame(cbind(time = Y_fit, status = delta_fit, Z_fit))
  colnames(df_cox) <- c("time", "status", z_names)

  formula_cox <- stats::as.formula(
    paste("survival::Surv(time, status) ~",
          paste(z_names, collapse = " + "))
  )

  cox_fit <- tryCatch(
    survival::coxph(formula_cox, data = df_cox, ties = "breslow"),
    error = function(e) NULL
  )

  if (is.null(cox_fit)) {
    beta_init        <- rep(0.0, ncol(Z))
    names(beta_init) <- z_names
  } else {
    beta_init        <- stats::coef(cox_fit)
    beta_init[is.na(beta_init)] <- 0.0
    names(beta_init) <- z_names
  }

  # At initialisation we assume all patients are uncured (W=1)
  W_unif  <- rep(1.0, n)
  breslow <- breslow_estimator(Y, delta, Z, beta_init, W_unif)
  Su_hat  <- pmax(eval_Su_diag(Y, Z, beta_init, breslow), 0.0)

  list(beta = beta_init, breslow = breslow, Su_hat = Su_hat)
}

initialize_em <- function(Y, delta, X, Z,
                          constraint = "gamma1",
                          verbose    = TRUE) {
  if (verbose) cat("  [1/2] Initialising incidence (logistic)...\n")
  inc <- init_incidence(X, delta, constraint)

  if (verbose) cat("  [2/2] Initialising latency (Cox on events only)...\n")
  lat <- init_latency(Y, delta, Z)

  if (verbose) {
    cat(sprintf("  Initial gamma : [%s]\n",
                paste(round(inc$gamma, 3), collapse=", ")))
    cat(sprintf("  Initial beta  : [%s]\n",
                paste(round(lat$beta,  3), collapse=", ")))
  }

  list(gamma   = inc$gamma,   p_hat     = inc$p_hat,
       beta    = lat$beta,    Su_hat    = lat$Su_hat,
       breslow = lat$breslow, h         = 0.7)
}
