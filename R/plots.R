# =============================================================
# plots.R
# Plotting functions for sic_fit objects
# =============================================================

#' Plot the Estimated Link Function of a Fitted SIC Model
#'
#' @description
#' Plots the estimated link function \eqn{\hat{g}(u)} against
#' the single index \eqn{\hat{\boldsymbol{\gamma}}^\top
#' \mathbf{x}_i} for the training observations.
#'
#' @param fit An object of class \code{"sic_fit"}.
#' @param main Character. Plot title. Default: automatic.
#' @param xlab Character. X-axis label. Default: automatic.
#' @param ylab Character. Y-axis label. Default: automatic.
#' @param col Color of the curve. Default: \code{"blue"}.
#' @param lwd Line width. Default: \code{2}.
#' @param ... Additional arguments passed to \code{plot}.
#'
#' @return Invisibly returns a data.frame with columns
#'   \code{index} and \code{g_hat}, sorted by index.
#'
#' @seealso \code{\link{sic_fit}}, \code{\link{plot_survival}}
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
#' n     <- 150
#' X     <- matrix(rnorm(n * 4), n, 4)
#' Z     <- matrix(rbinom(n, 1L, 0.6), n, 1)
#' p     <- 1 / (1 + exp(-X %*% c(1, -0.5, 0.3, -0.2)))
#' B     <- rbinom(n, 1L, p)
#' T0    <- ifelse(B == 1, rexp(n, 0.5), Inf)
#' C     <- rexp(n, 0.3)
#' Y     <- pmin(T0, C)
#' delta <- as.integer(T0 <= C & T0 < Inf)
#' fit   <- sic_fit(Y, delta, X, Z, h_step = 0.1, verbose = 0)
#' plot_link(fit)
#' }
#'
#' @importFrom graphics abline lines legend plot
#' @export
plot_link <- function(fit,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      col  = "blue",
                      lwd  = 2L,
                      ...) {

  if (!inherits(fit, "sic_fit"))
    stop("'fit' must be an object of class 'sic_fit'.")

  u_ord  <- order(fit$index_train)
  idx    <- fit$index_train[u_ord]
  g_hat  <- fit$p_hat[u_ord]

  if (is.null(main))
    main <- "Estimated link function g"
  if (is.null(xlab))
    xlab <- expression(hat(gamma)^T * x)
  if (is.null(ylab))
    ylab <- "Estimated uncure probability g(u)"

  plot(idx, g_hat,
       type = "l",
       col  = col,
       lwd  = lwd,
       ylim = c(0, 1),
       xlab = xlab,
       ylab = ylab,
       main = main,
       ...)

  abline(h = mean(fit$p_hat), lty = 2L,
         col = "gray60")

  invisible(data.frame(index = idx, g_hat = g_hat))
}


#' Plot the Estimated Latency Survival Function
#'
#' @description
#' Plots the estimated conditional survival function
#' \eqn{\hat{S}_u(t \mid \mathbf{z})} for the uncured subjects,
#' for specified values of the latency covariates \code{Z}.
#'
#' @param fit An object of class \code{"sic_fit"}.
#' @param newZ Numeric matrix of latency covariate values
#'   (\eqn{k \times q}), one row per group to plot. If
#'   \code{NULL} (default), uses \code{Z = 0} and \code{Z = 1}
#'   for a binary covariate.
#' @param t_eval Numeric vector of time points. Default:
#'   event times from training data.
#' @param group_labels Character vector of labels for each
#'   row of \code{newZ}. Default: \code{"Group 1"}, etc.
#' @param cols Color vector for each group. Default: automatic.
#' @param main Character. Plot title. Default: automatic.
#' @param xlab Character. X-axis label. Default: \code{"Time"}.
#' @param ylab Character. Y-axis label. Default: automatic.
#' @param lwd Line width. Default: \code{2}.
#' @param ... Additional arguments passed to \code{plot}.
#'
#' @return Invisibly returns a list with \code{t_eval} and
#'   \code{Su} matrix (one row per group).
#'
#' @seealso \code{\link{sic_fit}}, \code{\link{plot_link}}
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
#' n     <- 150
#' X     <- matrix(rnorm(n * 4), n, 4)
#' Z     <- matrix(rbinom(n, 1L, 0.6), n, 1)
#' p     <- 1 / (1 + exp(-X %*% c(1, -0.5, 0.3, -0.2)))
#' B     <- rbinom(n, 1L, p)
#' T0    <- ifelse(B == 1, rexp(n, 0.5), Inf)
#' C     <- rexp(n, 0.3)
#' Y     <- pmin(T0, C)
#' delta <- as.integer(T0 <= C & T0 < Inf)
#' fit   <- sic_fit(Y, delta, X, Z, h_step = 0.1, verbose = 0)
#'
#' # Plot for Z=0 and Z=1
#' plot_survival(fit,
#'               newZ = matrix(c(0, 1), ncol = 1),
#'               group_labels = c("Z = 0", "Z = 1"))
#' }
#'
#' @importFrom graphics lines legend plot
#' @export
plot_survival <- function(fit,
                          newZ         = NULL,
                          t_eval       = NULL,
                          group_labels = NULL,
                          cols         = NULL,
                          main         = NULL,
                          xlab         = "Time",
                          ylab         = NULL,
                          lwd          = 2L,
                          ...) {

  if (!inherits(fit, "sic_fit"))
    stop("'fit' must be an object of class 'sic_fit'.")

  # Default Z values : 0 and 1 for binary covariate
  if (is.null(newZ)) {
    q    <- ncol(fit$Z_train)
    newZ <- matrix(c(rep(0, q), rep(1, q)),
                   nrow = 2L, ncol = q, byrow = TRUE)
  }
  if (!is.matrix(newZ)) newZ <- as.matrix(newZ)

  k <- nrow(newZ)

  if (is.null(t_eval))
    t_eval <- fit$breslow$times
  if (length(t_eval) == 0L)
    stop("No event times available.")

  # Default labels
  if (is.null(group_labels))
    group_labels <- paste("Group", seq_len(k))

  # Default colors
  if (is.null(cols))
    cols <- c("blue","red","darkgreen","orange","purple")[
      seq_len(min(k, 5L))]
  if (length(cols) < k)
    cols <- rep(cols, length.out = k)

  # Compute Su for each group
  Su_mat <- eval_Su(t_eval, fit$breslow, newZ, fit$beta)

  if (is.null(ylab))
    ylab <- expression(hat(S)[u](t ~ "|" ~ z))
  if (is.null(main))
    main <- "Estimated latency survival S_u(t|z)"

  # Plot
  plot(t_eval, Su_mat[1L, ],
       type = "l", col = cols[1L], lwd = lwd,
       ylim = c(0, 1),
       xlab = xlab, ylab = ylab, main = main, ...)

  if (k > 1L) {
    for (i in seq(2L, k)) {
      lines(t_eval, Su_mat[i, ],
            col = cols[i], lwd = lwd)
    }
  }

  legend("topright",
         legend = group_labels,
         col    = cols[seq_len(k)],
         lwd    = lwd,
         bty    = "n")

  invisible(list(t_eval = t_eval, Su = Su_mat))
}
