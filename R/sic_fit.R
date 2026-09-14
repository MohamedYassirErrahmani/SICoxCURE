# =============================================================
# sic_fit.R
# Main user-facing function and S3 methods
# Reference: Amico, Van Keilegom & Legrand (2019)
#            Biometrics 75(2):452-462
# =============================================================

#' Fit the Single-Index/Cox (SIC) Mixture Cure Model
#'
#' @description
#' Fits the Single-Index/Cox (SIC) mixture cure model proposed
#' by Amico, Van Keilegom and Legrand (2019). The population
#' survival function is written as:
#'
#' \deqn{S(t \mid \mathbf{x}, \mathbf{z}) =
#'   1 - p(\mathbf{x}) + p(\mathbf{x}) \, S_u(t \mid \mathbf{z})}
#'
#' where \eqn{p(\mathbf{x}) = g(\boldsymbol{\gamma}^\top \mathbf{x})}
#' is the incidence (probability of being uncured), modeled via
#' estimated nonparametrically using kernel smoothing, and
#' \eqn{S_u(t \mid \mathbf{z})} is the latency (conditional
#' survival of the uncured subjects) modeled via a Cox
#' proportional hazards model.
#'
#' Estimation is performed by maximum likelihood via an
#' Expectation-Maximization (EM) algorithm. At each iteration,
#' the E-step computes the posterior probability of being
#' uncured (weights \eqn{W_i}), and the M-step maximizes the
#' expected complete-data log-likelihood separately for the
#' incidence and latency parts.
#'
#' @param Y Numeric vector of observed follow-up times
#'   \eqn{Y_i = \min(T_i, C_i)}, where \eqn{T_i} is the
#'   event time and \eqn{C_i} is the censoring time.
#'   Alternatively, a \code{data.frame} containing all
#'   variables (see \code{data} argument).
#'
#' @param delta Integer vector of event indicators
#'   (\code{1} = event observed, \code{0} = censored).
#'   Ignored if \code{data} is provided.
#'
#' @param X Numeric matrix of incidence covariates (\eqn{n
#'   \times d}). These are the covariates entering the
#'   single-index \eqn{\boldsymbol{\gamma}^\top \mathbf{X}}.
#'   Ignored if \code{data} is provided.
#'
#' @param Z Numeric matrix of latency covariates (\eqn{n
#'   \times q}). These enter the Cox model for the uncured
#'   subjects. Can overlap with \code{X}. If \code{NULL},
#'   \code{X} is used. Ignored if \code{data} is provided.
#'
#' @param data Optional \code{data.frame} containing all
#'   variables. If provided, \code{col_Y}, \code{col_delta},
#'   \code{cols_X} must be specified.
#'
#' @param col_Y Character. Name of the follow-up time column
#'   in \code{data}.
#'
#' @param col_delta Character. Name of the event indicator
#'   column in \code{data}.
#'
#' @param cols_X Character vector. Names of the incidence
#'   covariate columns in \code{data}.
#'
#' @param cols_Z Character vector. Names of the latency
#'   covariate columns in \code{data}. If \code{NULL},
#'   \code{cols_X} is used.
#'
#' @param kernel Character string specifying the kernel
#'   function for Nadaraya-Watson estimation of \eqn{g}.
#'   One of \code{"epanechnikov"} (default, as in Amico et
#'   al. 2019), \code{"gaussian"}, or \code{"uniform"}.
#'
#' @param h_range Numeric vector \code{c(h_min, h_max)}
#'   defining the search interval for bandwidth selection via
#'   likelihood cross-validation. Default: \code{c(0.4, 1.0)}
#'   as proposed by Amico, Van Keilegom and Legrand (2019) <doi:10.1111/biom.12999>.
#'
#' @param h_step Numeric. Step size of the bandwidth grid.
#'   Default: \code{0.001}, giving 601 candidate values in
#'   \code{[0.4, 1.0]}. Increase (e.g. \code{0.01}) for
#'   faster computation at some loss of precision.
#'
#' @param constraint Character string specifying the
#'   identifiability constraint on \eqn{\boldsymbol{\gamma}}.
#'   \code{"gamma1"} (default): \eqn{\gamma_1 = 1}.
#'   \code{"norm1"}: \eqn{\|\boldsymbol{\gamma}\| = 1} and
#'   \eqn{\gamma_1 > 0}. See Assumption A1-(v) in Amico et
#'   al. (2019).
#'
#' @param tol Numeric. Convergence tolerance for the EM
#'   algorithm. The algorithm stops when the maximum absolute
#'   change in \eqn{\boldsymbol{\gamma}}, \eqn{\boldsymbol{\beta}},
#'   and \eqn{S_0(\cdot)} between two consecutive iterations
#'   is smaller than \code{tol}. Default: \code{1e-5} as in
#'   Amico et al. (2019).
#'
#' @param max_iter Integer. Maximum number of EM iterations.
#'   Default: \code{500L}.
#'
#' @param conv_method Character string specifying how the
#'   convergence criterion is computed across parameters.
#'   \code{"both"} (default): minimum of the maximum and sum
#'   criteria. \code{"max"}: maximum absolute difference.
#'   \code{"sum"}: sum of absolute differences.
#'
#' @param standardize_X Logical. Should the incidence
#'   covariates \code{X} be standardized (zero mean, unit
#'   variance) before fitting? Default: \code{TRUE}.
#'   Standardization improves numerical stability and eases
#'   interpretation of \eqn{\boldsymbol{\gamma}} coefficients
#'   (Amico et al. 2019).
#'
#' @param verbose Integer controlling the verbosity level.
#'   \code{0}: silent. \code{1}: progress bar and convergence
#'   message (default). \code{2}: detailed output at each EM
#'   iteration (criterion, log-likelihood, bandwidth).
#'
#' @return An object of class \code{"sic_fit"}, which is a
#'   named list containing:
#'
#' \describe{
#'   \item{\code{gamma}}{Named numeric vector of length \eqn{d}.
#'     Estimated incidence coefficients
#'     \eqn{\hat{\boldsymbol{\gamma}}}. Under \code{constraint
#'     = "gamma1"}, \code{gamma[1] = 1} by construction.}
#'
#'   \item{\code{beta}}{Named numeric vector of length \eqn{q}.
#'     Estimated latency coefficients
#'     \eqn{\hat{\boldsymbol{\beta}}} from the Cox model.}
#'
#'   \item{\code{g_hat}}{Numeric vector of length \eqn{n}.
#'     Estimated uncure probabilities \eqn{\hat{p}(\mathbf{x}_i)}
#'     for the training observations.}
#'
#'   \item{\code{p_hat}}{Same as \code{g_hat} (alias for
#'     compatibility).}
#'
#'   \item{\code{Su_hat}}{Numeric vector of length \eqn{n}.
#'     Estimated latency survival \eqn{\hat{S}_u(Y_i \mid
#'     \mathbf{z}_i)} for each training observation at its
#'     own follow-up time.}
#'
#'   \item{\code{breslow}}{Named list with components
#'     \code{times} (ordered event times), \code{lambda0}
#'     (baseline hazard increments), \code{Lambda0}
#'     (cumulative baseline hazard), and \code{S0} (baseline
#'     survival function). Based on the weighted Breslow
#'     estimator (equation (5) of Amico et al. 2019).}
#'
#'   \item{\code{h}}{Numeric scalar. Selected bandwidth
#'     \eqn{\hat{h}} chosen by likelihood cross-validation
#'     (Amico et al. 2019).}
#'
#'   \item{\code{W}}{Numeric vector of length \eqn{n}. Final
#'     EM weights \eqn{W_i = \hat{P}(B_i = 1 \mid \text{data},
#'     \hat{\theta})}. Values in \eqn{[0, 1]}: \eqn{W_i = 1}
#'     for observed events, \eqn{W_i = 0} for censored
#'     observations beyond the cure threshold \eqn{\tau}.}
#'
#'   \item{\code{loglik}}{Numeric vector of length
#'     \code{n_iter}. Complete-data log-likelihood
#'     \eqn{\tilde{\ell}_c} at each EM iteration.}
#'
#'   \item{\code{loglik_comp}}{Same as \code{loglik} (alias).}
#'
#'   \item{\code{converged}}{Logical. \code{TRUE} if the EM
#'     algorithm converged within \code{max_iter} iterations.}
#'
#'   \item{\code{n_iter}}{Integer. Number of EM iterations
#'     performed.}
#'
#'   \item{\code{crit}}{Numeric. Final value of the
#'     convergence criterion.}
#'
#'   \item{\code{data_info}}{Named list with data summary:
#'     \code{n} (sample size), \code{d} (incidence dimension),
#'     \code{q} (latency dimension), \code{tau} (cure
#'     threshold), \code{n_events} (number of events),
#'     \code{n_plateau} (observations beyond \eqn{\tau}),
#'     \code{varnames_X}, \code{varnames_Z},
#'     \code{X_center}, \code{X_scale}.}
#'
#'   \item{\code{X_train}}{Numeric matrix. Standardized
#'     incidence covariates used for training.}
#'
#'   \item{\code{Z_train}}{Numeric matrix. Latency covariates
#'     used for training.}
#'
#'   \item{\code{Y_train}}{Numeric vector. Follow-up times
#'     used for training.}
#'
#'   \item{\code{delta_train}}{Integer vector. Event
#'     indicators used for training.}
#'
#'   \item{\code{index_train}}{Numeric vector of length
#'     \eqn{n}. Single indices
#'     \eqn{\hat{\boldsymbol{\gamma}}^\top \mathbf{x}_i} for
#'     the training observations.}
#'
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' The SIC cure model generalizes the classical logistic/Cox
#' for the incidence with an unknown function \eqn{g}
#' estimated nonparametrically. This avoids the S-shape
#' constraint of the logistic model while maintaining the
#' dimensionality reduction of the single-index structure.
#'
#' The bandwidth \eqn{h} is selected at each EM iteration by
#' minimizing the cross-validation criterion
#' \deqn{CV^{(m)}(h) = -\sum_{i=1}^n \left[
#'   W_i^{(m)} \log \hat{g}^{(m-1)}_{h,-i}(\hat{\boldsymbol{\gamma}}^{(m-1)\top}\mathbf{x}_i) +
#'   (1 - W_i^{(m)}) \log\left(1 - \hat{g}^{(m-1)}_{h,-i}(\hat{\boldsymbol{\gamma}}^{(m-1)\top}\mathbf{x}_i)\right)
#' \right]}
#' over a grid of \eqn{h} values in \code{h_range}, where
#' \eqn{\hat{g}_{h,-i}} is the leave-one-out Nadaraya-Watson
#' estimator with bandwidth \eqn{h}.
#'
#' Initial values are obtained from a logistic regression of
#' \code{delta} on \code{X} (for \eqn{\boldsymbol{\gamma}^{(0)}})
#' and a Cox model fitted on uncensored observations only
#' (for \eqn{\boldsymbol{\beta}^{(0)}}).
#'
#' @references
#' Amico M, Van Keilegom I, Legrand C (2019). The
#' single-index/Cox mixture cure model.
#' \emph{Biometrics}, \strong{75}(2), 452--462.
#' \doi{10.1111/biom.12999}
#'
#' Sy JP, Taylor JMG (2000). Estimation in a Cox proportional
#' hazards cure model. \emph{Biometrics}, \strong{56},
#' 227--236.
#'
#' Taylor JMG (1995). Semi-parametric estimation in failure
#' time mixture models. \emph{Biometrics}, \strong{51},
#' 899--907.
#'
#' Klein RW, Spady RH (1993). An efficient semiparametric
#' estimator for binary response models.
#' \emph{Econometrica}, \strong{61}, 387--421.
#'
#' @seealso
#' standard errors and p-values,
#' the SIC model,
#' prediction error.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n     <- 150
#' X     <- matrix(rnorm(n * 4), n, 4)
#' Z     <- matrix(rbinom(n, 1L, 0.6), n, 1)
#' p     <- 1 / (1 + exp(-X %*% c(1, -0.5, 0.3, -0.2)))
#' B     <- rbinom(n, 1L, p)
#' T0    <- ifelse(B == 1, rexp(n, 0.5), Inf)
#' C     <- rexp(n, 0.3)
#' Y     <- pmin(T0, C)
#' delta <- as.integer(T0 <= C & T0 < Inf)
#'
#' fit <- sic_fit(Y, delta, X, Z, h_step = 0.1, verbose = 0)
#' print(fit)
#' summary(fit)
#' }
#'
#' @export
sic_fit <- function(Y       = NULL,
                    delta   = NULL,
                    X       = NULL,
                    Z       = NULL,
                    data      = NULL,
                    col_Y     = NULL,
                    col_delta = NULL,
                    cols_X    = NULL,
                    cols_Z    = NULL,
                    kernel        = "epanechnikov",
                    h_range       = c(0.4, 1.0),
                    h_step        = 0.001,
                    constraint    = "gamma1",
                    tol           = 1e-5,
                    max_iter      = 500L,
                    conv_method   = c("both", "max", "sum"),
                    standardize_X = TRUE,
                    verbose       = 1L) {

  cl <- match.call()
  conv_method <- match.arg(conv_method)

  # -----------------------------------------------------------
  # data.frame interface: automatic column extraction
  # -----------------------------------------------------------
  if (!is.null(data)) {

    if (is.null(col_Y))
      stop("Argument 'col_Y' is missing: specify the name of the time column.")
    if (is.null(col_delta))
      stop("Argument 'col_delta' is missing: specify the name of the event column.")
    if (is.null(cols_X))
      stop("Argument 'cols_X' is missing: specify the names of the incidence covariate columns.")

    cols_missing <- setdiff(c(col_Y, col_delta, cols_X, cols_Z),
                            colnames(data))
    if (length(cols_missing) > 0L)
      stop(sprintf(
        "Columns not found in 'data': %s\n  Available columns: %s",
        paste(cols_missing, collapse = ", "),
        paste(colnames(data), collapse = ", ")))

    if (verbose >= 1L) {
      cat("=== SIC Model: extracting variables from data.frame ===\n")
      cat(sprintf("  Y     <- data$%s\n", col_Y))
      cat(sprintf("  delta <- data$%s\n", col_delta))
      cat(sprintf("  X     <- data[, c(%s)]\n",
                  paste(sprintf("'%s'", cols_X), collapse = ", ")))
      if (!is.null(cols_Z))
        cat(sprintf("  Z     <- data[, c(%s)]\n",
                    paste(sprintf("'%s'", cols_Z), collapse = ", ")))
      cat("\n")
    }

    Y     <- data[[col_Y]]
    delta <- data[[col_delta]]
    X     <- data[, cols_X, drop = FALSE]
    Z     <- if (!is.null(cols_Z)) data[, cols_Z, drop = FALSE] else NULL

  } else {

    if (is.null(Y))
      stop("Argument 'Y' is missing. Provide either Y/delta/X or data/col_Y/col_delta/cols_X.")
    if (is.null(delta))
      stop("Argument 'delta' is missing.")
    if (is.null(X))
      stop("Argument 'X' is missing.")
  }

  # -----------------------------------------------------------
  # Main pipeline
  # -----------------------------------------------------------
  if (verbose >= 1L) cat("=== Fitting SIC Mixture Cure Model ===\n")
  if (verbose >= 1L) cat("Step 1/3: Data validation and preparation...\n")

  prep <- prepare_data(Y, delta, X, Z, standardize_X)
  if (verbose >= 1L) print_data_summary(prep$data_info)

  if (verbose >= 1L) cat("Step 2/3: Initialisation...\n")
  init <- initialize_em(prep$Y, prep$delta, prep$X, prep$Z,
                        constraint, verbose = (verbose >= 1L))

  if (verbose >= 1L) cat("Step 3/3: EM algorithm...\n")
  em <- run_em(prep$Y, prep$delta, prep$X, prep$Z,
               init, prep$data_info,
               kernel, h_range, h_step,
               constraint, tol, max_iter,
               conv_method, verbose)

  names(em$gamma) <- prep$data_info$varnames_X
  names(em$beta)  <- prep$data_info$varnames_Z

  result <- c(em, list(
    call          = cl,
    kernel        = kernel,
    constraint    = constraint,
    h_range       = h_range,
    h_step        = h_step,
    tol           = tol,
    conv_method   = conv_method,
    data_info     = prep$data_info,
    X_train       = prep$X,
    Z_train       = prep$Z,
    Y_train       = prep$Y,
    delta_train   = prep$delta,
    index_train   = compute_index(prep$X, em$gamma),
    standardize_X = standardize_X
  ))

  class(result) <- "sic_fit"
  result
}

# =============================================================
# S3 methods
# =============================================================

#' Print a Fitted SIC Model
#'
#' @description
#' Prints a compact summary of a fitted \code{"sic_fit"} object,
#' showing the main model parameters and convergence status.
#'
#' @param x An object of class \code{"sic_fit"}.
#' @param digits Integer. Number of significant digits to print.
#'   Default: \code{4}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns \code{x}.
#'
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n     <- 150
#' X     <- matrix(rnorm(n * 4), n, 4)
#' Z     <- matrix(rbinom(n, 1L, 0.6), n, 1)
#' p     <- 1 / (1 + exp(-X %*% c(1, -0.5, 0.3, -0.2)))
#' B     <- rbinom(n, 1L, p)
#' T0    <- ifelse(B == 1, rexp(n, 0.5), Inf)
#' C     <- rexp(n, 0.3)
#' Y     <- pmin(T0, C)
#' delta <- as.integer(T0 <= C & T0 < Inf)
#'
#' fit <- sic_fit(Y, delta, X, Z, h_step = 0.1, verbose = 0)
#' print(fit)
#' summary(fit)
#' }
#'
#' @export
#' @method print sic_fit
print.sic_fit <- function(x, digits = 4L, ...) {
  cat("\n=== SIC Mixture Cure Model ===\n\n")
  cat(sprintf("Call: %s\n\n", deparse(x$call)))
  cat(sprintf("n = %d | events = %d | kernel = %s | h = %.4f\n",
              x$data_info$n,
              x$data_info$n_events,
              x$kernel,
              x$h))
  cat(sprintf("Convergence: %s in %d iterations (crit = %.2e)\n\n",
              ifelse(x$converged, "YES", "NO"),
              x$n_iter,
              x$crit))

  cat("Incidence coefficients (gamma):\n")
  print(round(x$gamma, digits), ...)

  cat("\nLatency coefficients (beta):\n")
  print(round(x$beta, digits), ...)

  cat("\nUse summary() for standard errors and p-values.\n")
  invisible(x)
}

#' Summarize a Fitted SIC Model
#'
#' @description
#' Produces a detailed summary of a fitted \code{"sic_fit"}
#' object, including coefficient estimates, standard errors,
#' z-statistics, and p-values (when standard errors are
#' provided via bootstrap).
#'
#' @param object An object of class \code{"sic_fit"}.
#' @param se Optional. Either a \code{"sic_bootstrap"} object
#'   returned by \code{\link{bootstrap_se}}, or a named list
#'   with components \code{gamma} and \code{beta} containing
#'   bootstrap standard errors. If \code{NULL} (default),
#'   only point estimates are shown.
#' @param digits Integer. Number of digits for rounding.
#'   Default: \code{4}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{object}.
#'
#' @details
#' Standard errors are not computed by default as they require
#' bootstrap resampling, which can be computationally intensive.
#' See the package vignette for an example of bootstrap
#' inference.
#'
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n     <- 150
#' X     <- matrix(rnorm(n * 4), n, 4)
#' Z     <- matrix(rbinom(n, 1L, 0.6), n, 1)
#' p     <- 1 / (1 + exp(-X %*% c(1, -0.5, 0.3, -0.2)))
#' B     <- rbinom(n, 1L, p)
#' T0    <- ifelse(B == 1, rexp(n, 0.5), Inf)
#' C     <- rexp(n, 0.3)
#' Y     <- pmin(T0, C)
#' delta <- as.integer(T0 <= C & T0 < Inf)
#'
#' fit <- sic_fit(Y, delta, X, Z, h_step = 0.1, verbose = 0)
#' print(fit)
#' summary(fit)
#' }
#'
#' @importFrom stats symnum
#' @export
#' @method summary sic_fit
summary.sic_fit <- function(object, se = NULL,
                            digits = 4L, ...) {
  di <- object$data_info

  cat("\n=== Summary: SIC Mixture Cure Model ===\n\n")
  cat(sprintf("Call: %s\n\n", deparse(object$call)))

  # Data summary
  cat("Data:\n")
  cat(sprintf("  n = %d | Events: %d (%.1f%%)",
              di$n, di$n_events, 100 * di$event_rate))
  cat(sprintf(" | Cure threshold tau = %.4f\n", di$tau))
  cat(sprintf("  Plateau: %d observations (%.1f%%)\n\n",
              di$n_plateau, di$pct_plateau))

  # Model parameters
  cat("Model:\n")
  cat(sprintf("  Kernel: %s | Bandwidth h = %.4f\n",
              object$kernel, object$h))
  cat(sprintf("  Constraint: '%s' | Standardized X: %s\n\n",
              object$constraint,
              ifelse(object$standardize_X, "yes", "no")))

  # Convergence
  cat(sprintf("Convergence: %s in %d iterations (crit = %.2e, tol = %.0e)\n\n",
              ifelse(object$converged, "YES", "NO"),
              object$n_iter, object$crit, object$tol))

  # Incidence table
  cat("Incidence (gamma) -- single-index coefficients:\n")
  tab_g <- data.frame(Estimate = round(object$gamma, digits),
                      row.names = names(object$gamma))
  if (!is.null(se) && !is.null(se$gamma)) {
    tab_g$Std.Error <- round(se$gamma, digits)
    tab_g$z.value   <- round(object$gamma / se$gamma, 3L)
    tab_g$Pr...z..  <- round(
      2 * stats::pnorm(-abs(object$gamma / se$gamma)), 4L)
    names(tab_g)[4L] <- "Pr(>|z|)"
    tab_g$sig <- symnum(tab_g$`Pr(>|z|)`,
                        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols   = c("***","**","*",".",""))
    names(tab_g)[5L] <- ""
  }
  print(tab_g)

  if (is.null(se))
    cat("(Standard errors not provided. Use bootstrap.)\n")

  # Latency table
  cat("\nLatency (beta) -- Cox model coefficients:\n")
  tab_b <- data.frame(Estimate = round(object$beta, digits),
                      row.names = names(object$beta))
  if (!is.null(se) && !is.null(se$beta)) {
    tab_b$Std.Error <- round(se$beta, digits)
    tab_b$z.value   <- round(object$beta / se$beta, 3L)
    tab_b$Pr...z..  <- round(
      2 * stats::pnorm(-abs(object$beta / se$beta)), 4L)
    names(tab_b)[4L] <- "Pr(>|z|)"
    tab_b$sig <- symnum(tab_b$`Pr(>|z|)`,
                        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols   = c("***","**","*",".",""))
    names(tab_b)[5L] <- ""
  }
  print(tab_b)

  if (is.null(se))
    cat("(Standard errors not provided. Use bootstrap.)\n")

  # Significance legend
  if (!is.null(se))
    cat("\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  invisible(object)
}
