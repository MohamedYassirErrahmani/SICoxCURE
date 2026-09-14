# =============================================================
# nadaraya_watson.R
# =============================================================

# ── Leave-one-out Nadaraya-Watson estimator ──────────────────
# training patient using a weighted average of W of neighbors
nw_loo <- function(u, W, h, kernel_fn) {
  n     <- length(u)
  diffs <- outer(u, u, "-") / h
  K_mat <- matrix(kernel_fn(as.vector(diffs)), n, n)

  # Leave-one-out : patient i cannot be her own neighbor
  diag(K_mat) <- 0.0

  num   <- as.vector(K_mat %*% W)

  denom <- rowSums(K_mat)

  # Detect patients with no nearby neighbor
  no_neighbors <- denom < 1e-10
  if (any(no_neighbors)) {
    warning(sprintf(
      "%d subject(s) have no nearby neighbor (h=%.3f may be too
       small, or covariates poorly standardized). Using mean(W)
       as fallback for these subjects.",
      sum(no_neighbors), h
    ))
  }

  g_hat <- ifelse(no_neighbors, mean(W), num / denom)
  clip_prob(g_hat)
}

# ── Nadaraya-Watson predictor for new patients ───────────────
# Predicts g for new patients using training patients as reference
# No leave-one-out needed : new patients ≠ training patients
nw_predict <- function(u_train, W_train, u_new, h, kernel_fn) {
  n_new <- length(u_new)
  n_tr  <- length(u_train)

  # RECTANGULAR matrix n_new x n_tr (vs SQUARE n x n in nw_loo)
  diffs <- outer(u_new, u_train, "-") / h
  K_mat <- matrix(kernel_fn(as.vector(diffs)), n_new, n_tr)

  num   <- as.vector(K_mat %*% W_train)
  denom <- rowSums(K_mat)

  no_neighbors <- denom < 1e-10
  if (any(no_neighbors)) {
    warning(sprintf(
      "%d new subject(s) have no nearby training neighbor
       (h=%.3f may be too small, or their covariate values
       lie outside the training range). Using mean(W_train)
       as fallback for these subjects.",
      sum(no_neighbors), h
    ))
  }

  g_hat <- ifelse(no_neighbors, mean(W_train), num / denom)
  clip_prob(g_hat)
}

# ── Cross-validation criterion for bandwidth selection ───────
# at gamma^(m-1)"
# = -log(L1_tilde) = quality score for ONE value of h
# SMALL score -> GOOD h, LARGE score -> BAD h
cv_criterion <- function(h, u, W, kernel_fn) {
  g_hat <- nw_loo(u, W, h, kernel_fn)
  val   <- -sum(W *safe_log(g_hat) +
                  (1.0 - W) *safe_log(1.0 - g_hat))
  if (is.nan(val) || is.infinite(val)) {
    warning(sprintf(
      "Cross-validation criterion is invalid (NaN or Infinite)
       for h=%.4f. Returning a large penalty value instead.",
      h
    ))
    return(1e15)
  }
  val
}

# ── Bandwidth selection via likelihood cross-validation ──────
# cross-validation criterion, computed at each EM iteration"
#
# Parameters :
#   u, W, kernel_fn : passed to cv_criterion()
#   h_range         : search interval [h_min, h_max]
#   h_step          : grid step size
#                     default 0.001 -> 601 values in [0.4, 1.0]
#                     user can set e.g. h_step=0.01 (61 values)
#                     for faster computation
select_bandwidth <- function(u, W, kernel_fn,
                             h_range = c(0.4, 1.0),
                             h_step  = 0.001) {

  # Build grid with step h_step
  h_grid <- seq(h_range[1], h_range[2], by = h_step)

  # Evaluate cv_criterion for each h in h_grid
  cv_vals <- vapply(h_grid,
                    function(h) cv_criterion(h, u, W, kernel_fn),
                    numeric(1L))

  # Find the h that minimises the CV criterion
  idx_opt <- which.min(cv_vals)

  list(h_opt   = h_grid[idx_opt],
       cv_vals = cv_vals,
       h_grid  = h_grid)
}
