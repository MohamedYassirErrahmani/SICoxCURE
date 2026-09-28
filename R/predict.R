# =============================================================
# predict.R
# Prediction method for a fitted SIC mixture cure model
# Reference: Amico, Van Keilegom & Legrand (2019)
#            Biometrics 75(2):452-462
# =============================================================

#' Predict from a Fitted SIC Mixture Cure Model
#'
#' @description
#' Computes predictions from a fitted \code{"sic_fit"} object.
#' Three types of predictions are available:
#' \itemize{
#'   \item \strong{Incidence}: estimated uncure probability
#'     \eqn{\hat{p}(\mathbf{x}_i) = \hat{g}(\hat{\boldsymbol{\gamma}}^\top
#'     \mathbf{x}_i)} for each observation.
#'   \item \strong{Latency}: estimated conditional survival
#'     \eqn{\hat{S}_u(t \mid \mathbf{z}_i)} for uncured subjects.
#'   \item \strong{Population survival}: estimated marginal survival
#'     \eqn{\hat{S}(t \mid \mathbf{x}_i, \mathbf{z}_i) =
#'     1 - \hat{p}(\mathbf{x}_i) + \hat{p}(\mathbf{x}_i) \,
#'     \hat{S}_u(t \mid \mathbf{z}_i)}.
#' }
#'
#' @param object An object of class \code{"sic_fit"}, as
#'   returned by \code{\link{sic_fit}}.
#'
#' @param newX Numeric matrix of new incidence covariates
#'   (\eqn{m \times d}). If \code{NULL} (default), predictions
#'   are returned for the training data.
#'
#' @param newZ Numeric matrix of new latency covariates
#'   (\eqn{m \times q}). If \code{NULL} and \code{newX} is
#'   provided, \code{newX} is used as \code{newZ}.
#'
#' @param t_eval Numeric vector of time points at which the
#'   latency and population survival are evaluated. Default:
#'   the observed event times from the training data.
#'
#' @param type Character string specifying the type of
#'   prediction. One of \code{"incidence"}, \code{"latency"},
#'   \code{"survival"}, or \code{"all"} (default).
#'
#' @param ... Additional arguments (currently ignored).
#'
#' @return Depends on \code{type}:
#' \describe{
#'   \item{\code{"incidence"}}{Numeric vector of length \eqn{m}.}
#'   \item{\code{"latency"}}{Matrix \eqn{m \times}
#'     \code{length(t_eval)}.}
#'   \item{\code{"survival"}}{Matrix \eqn{m \times}
#'     \code{length(t_eval)}.}
#'   \item{\code{"all"}}{Named list with \code{incidence},
#'     \code{latency}, \code{survival}, \code{t_eval},
#'     \code{index}.}
#' }
#'
#' @details
#' If \code{object$standardize_X} is \code{TRUE}, the columns
#' of \code{newX} are standardized using the mean and standard
#' deviation computed on the training data.
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
#'
#' fit <- sic_fit(Y, delta, X, Z, h_step = 0.1, verbose = 0)
#'
#' # Incidence on training data
#' p_hat <- predict(fit, type = "incidence")
#'
#' # All predictions on new data
#' pred <- predict(fit,
#'                 newX   = X[1:10, ],
#'                 newZ   = Z[1:10, , drop = FALSE],
#'                 type   = "all")
#' plot(pred$t_eval, pred$survival[1, ],
#'      type = "l", xlab = "Time", ylab = "S(t|x,z)")
#' }
#'
#' @export
#' @method predict sic_fit
predict.sic_fit <- function(object,
                            newX   = NULL,
                            newZ   = NULL,
                            t_eval = NULL,
                            type   = "all",
                            ...) {

  type <- match.arg(type,
                    c("incidence","latency","survival","all"))

  # ── Training data or new data ─────────────────────────────
  if (is.null(newX)) {
    X_new    <- object$X_train
    Z_new    <- object$Z_train
    on_train <- TRUE
  } else {
    if (!is.matrix(newX)) newX <- as.matrix(newX)
    if (object$standardize_X) {
      newX <- sweep(newX, 2L,
                    object$data_info$X_center, "-")
      newX <- sweep(newX, 2L,
                    object$data_info$X_scale,  "/")
    }
    X_new    <- newX
    Z_new    <- if (is.null(newZ)) X_new else {
      if (!is.matrix(newZ)) as.matrix(newZ) else newZ
    }
    on_train <- FALSE
  }

  if (is.null(t_eval))
    t_eval <- object$breslow$times
  if (length(t_eval) == 0L)
    stop("No event times available in 'object$breslow$times'.")

  kernel_fn <- get_kernel(object$kernel)

  # ── Incidence ─────────────────────────────────────────────
  if (type %in% c("incidence","survival","all")) {
    u_new <- compute_index(X_new, object$gamma)
    p_new <- if (on_train) {
      object$p_hat
    } else {
      nw_predict(object$index_train, object$W,
                 u_new, object$h, kernel_fn)
    }
  }

  # ── Latency ───────────────────────────────────────────────
  if (type %in% c("latency","survival","all")) {
    Su_new <- eval_Su(t_eval, object$breslow,
                      Z_new, object$beta)
  }

  # ── Population survival ───────────────────────────────────
  if (type %in% c("survival","all")) {
    S_new <- (1.0 - p_new) + p_new * Su_new
  }

  switch(type,
    incidence = p_new,
    latency   = Su_new,
    survival  = S_new,
    all = list(
      incidence = p_new,
      latency   = Su_new,
      survival  = S_new,
      t_eval    = t_eval,
      index     = compute_index(X_new, object$gamma)
    )
  )
}
