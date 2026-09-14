# =============================================================
# data_check.R
# Input data validation and preparation
# =============================================================

prepare_data <- function(Y, delta, X, Z = NULL,
                         standardize_X = FALSE) {
  #                 TRUE  -> center and scale each column of X
  #                          useful when X contains numeric variables
  #                          with large value ranges

  # -----------------------------------------------------------
  # 1. Convert factors/character columns to dummies
  # -----------------------------------------------------------
  .convert_to_matrix <- function(mat, name) {

    if (is.data.frame(mat)) {
      factor_cols <- sapply(mat, is.factor)
      char_cols   <- sapply(mat, is.character)
      need_dummy  <- factor_cols | char_cols

      if (any(need_dummy)) {
        col_names <- names(mat)[need_dummy]
        for (col in col_names) {
          vals    <- mat[[col]]
          if (is.character(vals)) vals <- as.factor(vals)
          levels_ <- levels(vals)
          n_lev   <- length(levels_)

          if (n_lev == 1L) {
            warning(sprintf(
              "Column '%s' in %s has only one level -- ignored.",
              col, name))
            mat[[col]] <- NULL
            next
          }

          if (n_lev == 2L) {
            mat[[col]] <- as.integer(vals == levels_[2L])
            message(sprintf(
              "%s: '%s' (binary) encoded as 0/1 -- reference='%s', active='%s'.",
              name, col, levels_[1L], levels_[2L]
            ))
          } else {
            ref   <- levels_[1L]
            dummy <- sapply(levels_[-1L], function(lv) as.integer(vals == lv))
            colnames(dummy) <- paste0(col, "_", levels_[-1L])
            mat[[col]] <- NULL
            mat <- cbind(mat, as.data.frame(dummy))
            message(sprintf(
              "%s: '%s' (%d levels) -> %d dummies, reference='%s'. Columns: %s",
              name, col, n_lev, n_lev - 1L, ref,
              paste(paste0(col, "_", levels_[-1L]), collapse=", ")
            ))
          }
        }
      }
      mat <- as.matrix(mat)
      storage.mode(mat) <- "double"
    } else {
      if (!is.matrix(mat)) mat <- as.matrix(mat)
      if (!is.numeric(mat))
        stop(sprintf(
          "%s must be numeric, a numeric data.frame or a numeric matrix.",
          name))
      storage.mode(mat) <- "double"
    }
    mat
  }

  # -----------------------------------------------------------
  # 2. Initial conversion
  # -----------------------------------------------------------
  X <- .convert_to_matrix(X, "X")
  if (!is.null(Z)) Z <- .convert_to_matrix(Z, "Z")

  # -----------------------------------------------------------
  # 3. Missing data handling
  # -----------------------------------------------------------
  n_total  <- length(Y)
  na_Y     <- is.na(Y)
  na_delta <- is.na(delta)
  na_X     <- apply(is.na(X), 1, any)
  na_Z     <- if (!is.null(Z)) apply(is.na(Z), 1, any) else rep(FALSE, n_total)
  na_any   <- na_Y | na_delta | na_X | na_Z
  n_exclus <- sum(na_any)

  if (n_exclus > 0L) {
    lignes_exclues <- which(na_any)
    message(sprintf(
      "MISSING DATA: %d subject(s) excluded out of %d (%.1f%%).",
      n_exclus, n_total, 100.0 * n_exclus / n_total
    ))
    message(sprintf("  Excluded rows: %s",
      if (length(lignes_exclues) <= 10)
        paste(lignes_exclues, collapse=", ")
      else
        paste0(paste(lignes_exclues[1:10], collapse=", "),
               " ... (", length(lignes_exclues), " total)")
    ))
    ok    <- !na_any
    Y     <- Y[ok]
    delta <- delta[ok]
    X     <- X[ok, , drop=FALSE]
    if (!is.null(Z)) Z <- Z[ok, , drop=FALSE]
    message(sprintf("  => Analysis conducted on %d complete subjects.\n", sum(ok)))
  }

  n <- length(Y)

  # -----------------------------------------------------------
  # 4. Post-cleaning validations
  # -----------------------------------------------------------
  if (!is.numeric(Y) || any(Y < 0))
    stop("Y must be numeric and >= 0.")
  if (any(!delta %in% c(0L, 1L)))
    stop("delta must contain only 0 or 1.")
  if (length(delta) != n)
    stop("Y and delta must have the same length.")
  if (n < 20L)
    stop(sprintf("Too few complete subjects (%d). Minimum required: 20.", n))
  if (nrow(X) != n)
    stop("X must have as many rows as Y.")

  # -----------------------------------------------------------
  # 5. Column names and X validation
  # -----------------------------------------------------------
  varnames_X <- colnames(X)
  if (is.null(varnames_X)) {
    varnames_X  <- paste0("X", seq_len(ncol(X)))
    colnames(X) <- varnames_X
  }
  d <- ncol(X)

  X_vars <- apply(X, 2, stats::var)
  if (any(X_vars < 1e-10))
    stop(paste("Column(s) of X with zero variance:",
               paste(varnames_X[X_vars < 1e-10], collapse=", ")))

  # -----------------------------------------------------------
  # 6. Warning if numeric variables are not standardized
  #    Epanechnikov kernel works correctly.
  #    K(u) = 0 if u^2 > 5 -> if u_i values are large,
  #    all distances will be > sqrt(5) -> K=0 everywhere
  # -----------------------------------------------------------
  X_sds <- apply(X, 2, stats::sd)

  if (!standardize_X) {
    for (j in seq_len(d)) {
      sd_j    <- X_sds[j]
      range_j <- diff(range(X[, j]))

      # Warn if standard deviation is large (unstandardized numeric variable)
      # Threshold = 2 because after standardization sd = 1
      # If sd > 2 -> variable is probably not standardized
      if (sd_j > 2.0) {
        warning(sprintf(
          "X: column '%s' has standard deviation %.2f and range %.2f.
           The Epanechnikov kernel has support sqrt(5)=2.236.
           With large u_i values, K(u)=0 everywhere ->
           nw_loo() returns mean(W) for all subjects -> g poorly estimated.
           Article recommends standardization.
           Solution 1 (recommended): standardize manually before sic_fit()
              dat$%s_std <- scale(dat$%s)
           Solution 2: use standardize_X=TRUE in sic_fit()",
          varnames_X[j], sd_j, range_j,
          varnames_X[j], varnames_X[j]
        ))
      }
    }
  }

  # -----------------------------------------------------------
  # 7. Standardize X if requested
  #    continuous and discrete, is considered in the incidence"
  # -----------------------------------------------------------
  X_center <- NULL  # NULL if no standardization
  X_scale  <- NULL  # NULL if no standardization

  if (standardize_X) {
    X_center <- colMeans(X)
    X_scale  <- X_sds

    # Protection: if sd = 0 for a column (constant variable)
    # -> already detected in step 5 -> no need to handle here

    X <- scale(X, center = X_center, scale = X_scale)

    message(sprintf(
      "X standardized (centered and scaled): %d column(s). [Article]",
      d
    ))

    for (j in seq_len(d)) {
      message(sprintf(
        "  %s: mean=%.3f, sd=%.3f -> after: mean=0, sd=1",
        varnames_X[j], X_center[j], X_scale[j]
      ))
    }
  }

  # -----------------------------------------------------------
  # 8. Z preparation
  # -----------------------------------------------------------
  if (is.null(Z)) {
    Z <- X
    message("Z not provided: using same covariates as X for the latency.")
  }
  if (nrow(Z) != n)
    stop("Z must have as many rows as Y.")

  varnames_Z <- colnames(Z)
  if (is.null(varnames_Z)) {
    varnames_Z  <- paste0("Z", seq_len(ncol(Z)))
    colnames(Z) <- varnames_Z
  }
  q <- ncol(Z)

  # -----------------------------------------------------------
  # 9. Cure threshold
  # -----------------------------------------------------------
  tau        <- compute_cure_threshold(Y, delta)
  beyond_tau <- is_beyond_threshold(Y, delta, tau)
  n_plateau  <- sum(beyond_tau)
  n_events   <- sum(delta == 1L)

  if (n_plateau == 0L)
    warning("No subject beyond cure threshold tau = ", round(tau, 4))
  if (n_plateau < 5L)
    warning("Only ", n_plateau, " subject(s) in the plateau.")

  # -----------------------------------------------------------
  # 10. Build data_info
  # -----------------------------------------------------------
  data_info <- list(
    n            = n,
    n_initial    = n_total,
    n_exclus     = n_exclus,
    n_events     = n_events,
    n_censored   = n - n_events,
    event_rate   = n_events / n,
    tau          = tau,
    n_plateau    = n_plateau,
    pct_plateau  = 100.0 * n_plateau / n,
    d            = d,
    q            = q,
    varnames_X   = varnames_X,
    varnames_Z   = varnames_Z,
    beyond_tau   = beyond_tau,

    # Standardization info
    # NULL if standardize_X = FALSE
    # Used in predict.R to reverse standardization
    standardize_X = standardize_X,
    X_center      = X_center,
    X_scale       = X_scale
  )

  list(Y=Y, delta=as.integer(delta), X=X, Z=Z, data_info=data_info)
}

# -----------------------------------------------------------
# Data summary display
# -----------------------------------------------------------
print_data_summary <- function(data_info) {
  cat("=== Data summary ===\n")
  if (data_info$n_exclus > 0L) {
    cat(sprintf("  Initial sample: %d | Excluded (NA): %d | Analysed: %d\n",
                data_info$n_initial, data_info$n_exclus, data_info$n))
  } else {
    cat(sprintf("  n=%d (no missing data)\n", data_info$n))
  }
  cat(sprintf("  Events=%d (%.1f%%) | Censored=%d\n",
              data_info$n_events, 100.0 * data_info$event_rate, data_info$n_censored))
  cat(sprintf("  tau=%.4f | Plateau=%d (%.1f%%)\n",
              data_info$tau, data_info$n_plateau, data_info$pct_plateau))
  cat(sprintf("  X [%s] | Z [%s]\n",
              paste(data_info$varnames_X, collapse=", "),
              paste(data_info$varnames_Z, collapse=", ")))

  if (isTRUE(data_info$standardize_X)) {
    cat(sprintf("  X standardization: YES [Article]\n"))
  } else {
    cat(sprintf("  X standardization: NO (default)\n"))
  }

  cat("====================\n\n")
}
