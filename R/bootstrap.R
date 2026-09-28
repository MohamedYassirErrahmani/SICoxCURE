# =============================================================
# bootstrap.R
# Bootstrap standard errors for the SICoxCURE model
# Reference: Amico, Van Keilegom & Legrand (2019)
#            Biometrics 75(2):452-462
# =============================================================

#' Bootstrap Standard Errors for the SIC Mixture Cure Model
#'
#' @description
#' Computes bootstrap standard errors for the incidence
#' coefficients \eqn{\hat{\boldsymbol{\gamma}}} and the
#' latency coefficients \eqn{\hat{\boldsymbol{\beta}}} of a
#' fitted SIC mixture cure model, as in Amico, Van Keilegom
#' and Legrand (2019) \doi{10.1111/biom.12999}.
#'
#' At each bootstrap iteration, a sample of size \eqn{n} is
#' drawn with replacement from the training data, the SIC
#' model is re-fitted on this sample, and the resulting
#' estimates are stored. The standard errors are computed as
#' the standard deviation of the bootstrap estimates across
#' all \code{B} iterations.
#'
#' @param fit An object of class \code{"sic_fit"} returned by
#'   \code{\link{sic_fit}}.
#'
#' @param B Integer. Number of bootstrap samples. Default:
#'   \code{250} as in Amico, Van Keilegom and Legrand (2019)
#'   \doi{10.1111/biom.12999}.
#'
#' @param parallel Logical. Should the bootstrap be run in
#'   parallel using multiple CPU cores? Default: \code{TRUE}.
#'   Requires the \code{parallel} package (included in base R).
#'
#' @param ncores Integer. Number of CPU cores to use when
#'   \code{parallel = TRUE}. Default: \code{NULL}, which
#'   uses all available cores minus one
#'   (\code{parallel::detectCores() - 1}).
#'
#' @param seed Integer or \code{NULL}. Random seed for
#'   reproducibility. Default: \code{NULL}.
#'
#' @param verbose Logical. Should a progress message be
#'   printed? Default: \code{TRUE}.
#'
#' @return A named list with class \code{"sic_bootstrap"}
#'   containing:
#' \describe{
#'   \item{\code{gamma}}{Named numeric vector of length
#'     \eqn{d}. Bootstrap standard errors for the incidence
#'     coefficients \eqn{\hat{\boldsymbol{\gamma}}}.}
#'   \item{\code{beta}}{Named numeric vector of length
#'     \eqn{q}. Bootstrap standard errors for the latency
#'     coefficients \eqn{\hat{\boldsymbol{\beta}}}.}
#'   \item{\code{gamma_boot}}{Numeric matrix of dimensions
#'     \eqn{B \times d}. All bootstrap estimates of
#'     \eqn{\boldsymbol{\gamma}}.}
#'   \item{\code{beta_boot}}{Numeric matrix of dimensions
#'     \eqn{B \times q}. All bootstrap estimates of
#'     \eqn{\boldsymbol{\beta}}.}
#'   \item{\code{B}}{Integer. Number of bootstrap samples
#'     actually used (excluding failed iterations).}
#'   \item{\code{B_requested}}{Integer. Number of bootstrap
#'     samples requested.}
#'   \item{\code{n_failed}}{Integer. Number of bootstrap
#'     iterations that failed (e.g. convergence issues).}
#' }
#'
#' @details
#' The standard errors are then passed to
#' \code{\link{summary.sic_fit}} via the \code{se} argument
#' to obtain a coefficient table with p-values, as in Table 3
#' of Amico, Van Keilegom and Legrand (2019)
#' \doi{10.1111/biom.12999}.
#'
#' When \code{parallel = TRUE}, the function uses
#' \code{parallel::mclapply} on Unix/macOS and
#' \code{parallel::parLapply} on Windows. On Windows,
#' parallelization may be slower due to overhead; consider
#' setting \code{parallel = FALSE} for small datasets.
#'
#' Failed bootstrap iterations (non-convergence, numerical
#' issues) are automatically excluded. A warning is issued
#' if more than 10\% of iterations fail.
#'
#' @seealso
#' \code{\link{sic_fit}} for fitting the model,
#' \code{\link{summary.sic_fit}} for the coefficient table
#' with standard errors and p-values.
#'
#' @references
#' Amico M, Van Keilegom I, Legrand C (2019).
#' The single-index/Cox mixture cure model.
#' \emph{Biometrics}, \strong{75}(2), 452--462.
#' \doi{10.1111/biom.12999}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n     <- 200
#' X     <- matrix(rnorm(n * 4), n, 4)
#' Z     <- matrix(rbinom(n, 1L, 0.6), n, 1)
#' p     <- 1 / (1 + exp(-X %*% c(1, -0.5, 0.3, -0.2)))
#' B_ind <- rbinom(n, 1L, p)
#' T0    <- ifelse(B_ind == 1, rexp(n, 0.5), Inf)
#' C     <- rexp(n, 0.3)
#' Y     <- pmin(T0, C)
#' delta <- as.integer(T0 <= C & T0 < Inf)
#'
#' # Fit the model
#' fit <- sic_fit(Y, delta, X, Z, verbose = 0)
#'
#' # Compute bootstrap standard errors (B=50 for speed)
#' se <- bootstrap_se(fit, B = 50, parallel = FALSE,
#'                    seed = 1, verbose = TRUE)
#'
#' # Summary with p-values (as in Table 3 of the article)
#' summary(fit, se = se)
#'
#' # Access bootstrap distributions
#' hist(se$gamma_boot[, 2],
#'      main = "Bootstrap distribution of gamma[2]",
#'      xlab = "gamma[2]")
#' }
#'
#' @export
bootstrap_se <- function(fit,
                         B        = 250L,
                         parallel = TRUE,
                         ncores   = NULL,
                         seed     = NULL,
                         verbose  = TRUE) {

  if (!inherits(fit, "sic_fit"))
    stop("'fit' must be an object of class 'sic_fit'.")
  if (B < 2L)
    stop("'B' must be at least 2.")

  if (!is.null(seed)) set.seed(seed)

  # ── Data from training ───────────────────────────────────────
  Y      <- fit$Y_train
  delta  <- fit$delta_train
  X      <- fit$X_train
  Z      <- fit$Z_train
  n      <- length(Y)

  d <- length(fit$gamma)
  q <- length(fit$beta)

  if (verbose) {
    cat(sprintf(
      "Bootstrap standard errors: B=%d, n=%d, parallel=%s\n",
      B, n, ifelse(parallel, "TRUE", "FALSE")
    ))
    if (parallel) {
      nc <- if (is.null(ncores))
        max(1L, parallel::detectCores() - 1L)
      else as.integer(ncores)
      cat(sprintf("Using %d cores\n", nc))
    }
    cat("This may take several minutes...\n")
  }

  # ── Single bootstrap iteration ────────────────────────────────
  one_boot <- function(b) {
    idx <- sample(n, n, replace = TRUE)
    tryCatch({
      fit_b <- sic_fit(
        Y       = Y[idx],
        delta   = delta[idx],
        X       = X[idx, , drop = FALSE],
        Z       = Z[idx, , drop = FALSE],
        kernel        = fit$kernel,
        h_range       = fit$h_range,
        h_step        = fit$h_step,
        constraint    = fit$constraint,
        tol           = fit$tol,
        max_iter      = fit$data_info$max_iter %||% 500L,
        conv_method   = fit$conv_method,
        standardize_X = FALSE,
        verbose       = 0L
      )
      list(gamma = fit_b$gamma, beta = fit_b$beta,
           ok = TRUE)
    }, error = function(e) {
      list(gamma = rep(NA_real_, d),
           beta  = rep(NA_real_, q),
           ok    = FALSE)
    })
  }

  # ── Run bootstrap ─────────────────────────────────────────────
  seeds_b <- sample.int(.Machine$integer.max, B)

  if (parallel && .Platform$OS.type != "windows") {
    nc      <- if (is.null(ncores))
      max(1L, parallel::detectCores() - 1L)
    else as.integer(ncores)
    results <- parallel::mclapply(
      seq_len(B),
      function(b) { set.seed(seeds_b[b]); one_boot(b) },
      mc.cores = nc
    )
  } else if (parallel && .Platform$OS.type == "windows") {
    nc   <- if (is.null(ncores))
      max(1L, parallel::detectCores() - 1L)
    else as.integer(ncores)
    cl   <- parallel::makeCluster(nc)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("Y","delta","X","Z","fit","d","q",
                  "seeds_b","sic_fit","one_boot"),
      envir   = environment()
    )
    parallel::clusterEvalQ(cl, library(SICoxCURE))
    results <- parallel::parLapply(
      cl, seq_len(B),
      function(b) { set.seed(seeds_b[b]); one_boot(b) }
    )
  } else {
    # Sequential
    results <- vector("list", B)
    if (verbose) {
      pb <- utils::txtProgressBar(min=0, max=B, style=3)
    }
    for (b in seq_len(B)) {
      set.seed(seeds_b[b])
      results[[b]] <- one_boot(b)
      if (verbose) utils::setTxtProgressBar(pb, b)
    }
    if (verbose) close(pb)
  }

  # ── Collect results ───────────────────────────────────────────
  ok_idx <- which(sapply(results, function(r) r$ok))
  n_fail <- B - length(ok_idx)

  if (n_fail > 0) {
    warning(sprintf(
      "%d out of %d bootstrap iterations failed and were excluded.",
      n_fail, B
    ))
    if (n_fail > 0.1 * B)
      warning("More than 10%% of bootstrap iterations failed. Results may be unreliable.")
  }

  if (length(ok_idx) < 2L)
    stop("Fewer than 2 bootstrap iterations succeeded. Cannot compute standard errors.")

  gamma_boot <- do.call(rbind, lapply(results[ok_idx],
                                       `[[`, "gamma"))
  beta_boot  <- do.call(rbind, lapply(results[ok_idx],
                                       `[[`, "beta"))

  colnames(gamma_boot) <- names(fit$gamma)
  colnames(beta_boot)  <- names(fit$beta)

  # ── Standard errors = SD of bootstrap estimates ───────────────
  se_gamma <- apply(gamma_boot, 2L, stats::sd)
  se_beta  <- apply(beta_boot,  2L, stats::sd)

  if (verbose) {
    cat(sprintf(
      "\nDone. %d/%d iterations succeeded.\n",
      length(ok_idx), B
    ))
  }

  structure(
    list(
      gamma       = se_gamma,
      beta        = se_beta,
      gamma_boot  = gamma_boot,
      beta_boot   = beta_boot,
      B           = length(ok_idx),
      B_requested = B,
      n_failed    = n_fail
    ),
    class = "sic_bootstrap"
  )
}

# ── Null-coalescing operator ───────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Print Bootstrap Standard Errors
#'
#' @param x An object of class \code{"sic_bootstrap"}.
#' @param digits Integer. Number of digits. Default: \code{4}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.sic_bootstrap <- function(x, digits = 4L, ...) {
  cat("\n=== Bootstrap Standard Errors (SICoxCURE) ===\n\n")
  cat(sprintf("Bootstrap samples : %d / %d succeeded\n",
              x$B, x$B_requested))
  if (x$n_failed > 0)
    cat(sprintf("Failed iterations : %d\n", x$n_failed))
  cat("\nIncidence (gamma):\n")
  print(round(x$gamma, digits))
  cat("\nLatency (beta):\n")
  print(round(x$beta, digits))
  cat("\nUse summary(fit, se = <this object>) for p-values.\n")
  invisible(x)
}
