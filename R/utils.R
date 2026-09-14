# =============================================================
# utils.R
# Fonctions de base partagées par tout le package
# =============================================================

# --- NOYAUX ---

kernel_epanechnikov <- function(u) {
  val <- (3.0 / (4.0 * sqrt(5.0))) * (1.0 - u^2 / 5.0)
  val[u^2 > 5.0] <- 0.0
  pmax(val, 0.0)
}

kernel_gaussian <- function(u) {
  stats::dnorm(u)
}

kernel_uniform <- function(u) {
  0.5 * (abs(u) <= 1.0)
}

get_kernel <- function(kernel = c("epanechnikov", "gaussian", "uniform")) {
  kernel <- match.arg(kernel)
  switch(kernel,
         epanechnikov = kernel_epanechnikov,
         gaussian     = kernel_gaussian,
         uniform      = kernel_uniform
  )
}

# --- STABILITE NUMERIQUE ---

clip_prob <- function(p, eps = 1e-10) {
  pmax(eps, pmin(1.0 - eps, p))
}

safe_log <- function(x, eps = 1e-10) {
  log(pmax(x, eps))
}

safe_exp <- function(x, max_val = 700.0) {
  exp(pmin(x, max_val))
}

# --- INDICE UNIQUE ---

compute_index <- function(X, gamma) {
  if (!is.matrix(X)) X <- as.matrix(X)
  as.vector(X %*% gamma)
}

normalize_gamma <- function(gamma,
                            constraint = c("gamma1", "norm1")) {
  constraint <- match.arg(constraint)
  if (constraint == "gamma1") {
    if (abs(gamma[1]) < 1e-12) {
      warning("gamma[1] proche de 0.")
      return(list(gamma = gamma, scale = 1.0))
    }
    scale <- gamma[1]
    gamma_norm <- gamma / scale
  } else {
    nrm <- sqrt(sum(gamma^2))
    if (nrm < 1e-12) {
      return(list(gamma = gamma, scale = 1.0))
    }
    gamma_norm <- gamma / nrm
    if (gamma_norm[1] < 0.0) gamma_norm <- -gamma_norm
    scale <- nrm
  }
  list(gamma = gamma_norm, scale = scale)
}

# --- SEUIL DE GUERISON ---

compute_cure_threshold <- function(Y, delta) {
  event_times <- Y[delta == 1]
  if (length(event_times) == 0L)
    stop("Aucun evenement observe.")
  max(event_times)
}

is_beyond_threshold <- function(Y, delta, tau = NULL) {
  if (is.null(tau)) tau <- compute_cure_threshold(Y, delta)
  (delta == 0L) & (Y > tau)
}

# --- ESTIMATEUR DE BRESLOW PONDERE ---

breslow_estimator <- function(Y, delta, Z, beta, W) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)

  event_times <- sort(unique(Y[delta == 1L]))
  r           <- length(event_times)

  if (r == 0L) {
    return(list(times   = numeric(0), lambda0 = numeric(0),
                Lambda0 = numeric(0), S0      = numeric(0)))
  }

  lin_pred   <- as.vector(Z %*% beta)
  lin_pred   <- lin_pred - max(lin_pred)
  risk_score <- W *safe_exp(lin_pred)

  lambda0 <- numeric(r)
  for (j in seq_len(r)) {
    t_j        <- event_times[j]
    D_j        <- sum(delta[Y == t_j])
    denom      <- sum(risk_score[Y >= t_j])
    lambda0[j] <- if (denom > 1e-15) D_j / denom else 0.0
  }

  Lambda0 <- cumsum(lambda0)
  S0      <- exp(-Lambda0)

  list(times = event_times, lambda0 = lambda0,
       Lambda0 = Lambda0, S0 = S0)
}

eval_Su_diag <- function(Y, Z, beta, breslow) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  if (length(breslow$times) == 0L) return(rep(0.0, length(Y)))

  Lambda0_Yi <- stats::approx(
    x = c(0.0, breslow$times), y = c(0.0, breslow$Lambda0),
    xout = Y, method = "constant", rule = 2L, f = 0.0
  )$y

  lin_pred <- as.vector(Z %*% beta)
  lin_pred <- lin_pred - max(lin_pred)
  exp(-Lambda0_Yi *safe_exp(lin_pred))
}

eval_Su <- function(t_eval, breslow, Z_new, beta) {
  if (!is.matrix(Z_new)) Z_new <- as.matrix(Z_new)
  if (length(breslow$times) == 0L)
    return(matrix(0.0, nrow(Z_new), length(t_eval)))

  Lambda0_eval <- stats::approx(
    x = c(0.0, breslow$times), y = c(0.0, breslow$Lambda0),
    xout = t_eval, method = "constant", rule = 2L, f = 0.0
  )$y

  lin_pred <- as.vector(Z_new %*% beta)
  lin_pred <- lin_pred - max(lin_pred)
  outer(safe_exp(lin_pred), Lambda0_eval, function(e, L) exp(-e * L))
}

# --- CONVERGENCE ---
# of beta, gamma, and S0(.) is smaller than 10^-5"

convergence_criterion <- function(beta_new, beta_old,
                                  gamma_new, gamma_old,
                                  S0_new, S0_old,
                                  method = c("both", "max", "sum")) {
  method <- match.arg(method)

  # Differences absolues
  diff_beta  <- abs(beta_new  - beta_old)
  diff_gamma <- abs(gamma_new - gamma_old)
  diff_S0    <- if (length(S0_new) == length(S0_old))
    abs(S0_new - S0_old)
  else rep(Inf, 1L)

  # Critere MAX
  crit_max <- max(max(diff_beta),
                  max(diff_gamma),
                  max(diff_S0))

  # Critere SOMME
  crit_sum <- sum(sum(diff_beta),
                  sum(diff_gamma),
                  sum(diff_S0))

  switch(method,
         "both" = {
           crit <- min(crit_max, crit_sum)
           list(
             crit        = crit,
             crit_max    = crit_max,
             crit_sum    = crit_sum,
             method_used = if (crit_max < crit_sum) "max" else "sum"
           )
         },
         "max" = list(
           crit        = crit_max,
           crit_max    = crit_max,
           crit_sum    = crit_sum,
           method_used = "max"
         ),
         "sum" = list(
           crit        = crit_sum,
           crit_max    = crit_max,
           crit_sum    = crit_sum,
           method_used = "sum"
         )
  )
}

# =============================================================
# Vraisemblance COMPLETE L_c tilde UNIQUEMENT
# L_c = produit_i {p*lambda_u*Su}^(Wi*delta_i)
#               * {p*Su}^(Wi*(1-delta_i))
#               * {1-p}^((1-Wi)*(1-delta_i))
# =============================================================

complete_loglik <- function(p_hat, Su_hat, hu_hat, delta, W) {

  p_hat  <- clip_prob(p_hat)
  Su_hat <- pmax(Su_hat, 0.0)
  hu_hat <- pmax(hu_hat, 1e-15)
  W      <- pmax(0.0, pmin(1.0, W))

  # Part 1: subjects with observed event (delta=1)
  # {p(Xi)*lambda_u(Yi|Zi)*Su(Yi|Zi)}^(Wi*delta_i)
  ll_ev <- sum(W * delta *
               (safe_log(p_hat) +
                safe_log(hu_hat) +
                safe_log(Su_hat)))

  # Part 2: censored subjects (delta=0)
  # [{p(Xi)*Su(Yi|Zi)}^Wi * {1-p(Xi)}^(1-Wi)]^(1-delta_i)
  ll_cens <- sum((1L - delta) *
                 (W       * (safe_log(p_hat) + safe_log(Su_hat)) +
                  (1L - W) * safe_log(1.0 - p_hat)))

  ll_ev + ll_cens
}

# --- HAZARD DE LA LATENCE ---
# Required for complete_loglik()

compute_hazard_u <- function(Y, delta, Z, beta, breslow) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  n  <- length(Y)
  hu <- numeric(n)

  if (length(breslow$times) == 0L) return(pmax(hu, 1e-15))

  lin_pred <- as.vector(Z %*% beta) - max(as.vector(Z %*% beta))

  for (i in which(delta == 1L)) {
    j <- which(abs(breslow$times - Y[i]) < 1e-10)

    if (length(j) == 0L)
      stop(sprintf(
        "Patient %d: event time %.4f not found in breslow$times.",
        i, Y[i]))

    hu[i] <- breslow$lambda0[j[1L]] *safe_exp(lin_pred[i])

    if (hu[i] < 1e-15 && breslow$lambda0[j[1L]] >= 1e-15)
      warning(sprintf(
        "Patient %d: lambda0=%.6f but hu=0 (lin_pred=%.2f).",
        i, breslow$lambda0[j[1L]], lin_pred[i]))
  }

  pmax(hu, 1e-15)
}
