# =============================================================
# predict.R
# Prediction methods for a fitted SIC mixture cure model
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
#'     \mathbf{x}_i)} for each new observation.
#'   \item \strong{Latency}: estimated conditional survival
#'     \eqn{\hat{S}_u(t \mid \mathbf{z}_i)} for the uncured subjects.
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
#'   \item{\code{"latency"}}{Matrix \eqn{m \times} \code{length(t_eval)}.}
#'   \item{\code{"survival"}}{Matrix \eqn{m \times} \code{length(t_eval)}.}
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
#' @seealso \code{\link{sic_fit}}, \code{\link{prediction_error}}
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
#' ## Predictions on training data
#' p_hat <- predict(fit, type = "incidence")
#'
#' ## Predictions on new data
#' pred <- predict(fit,
#'                 newX = X[1:10, ],
#'                 newZ = Z[1:10, , drop = FALSE],
#'                 type = "all")
#' plot(pred$t_eval, pred$survival[1, ], type = "l",
#'      xlab = "Time", ylab = "S(t|x,z)")
#' }
#'
#' @export
#' @method predict sic_fit
predict.sic_fit <- function(object, newX = NULL, newZ = NULL,
                            t_eval = NULL, type = "all", ...) {

  type <- match.arg(type, c("incidence", "latency",
                             "survival", "all"))

  if (is.null(newX)) {
    X_new    <- object$X_train
    Z_new    <- object$Z_train
    on_train <- TRUE
  } else {
    if (!is.matrix(newX)) newX <- as.matrix(newX)
    if (object$standardize_X) {
      newX <- sweep(newX, 2L, object$data_info$X_center, "-")
      newX <- sweep(newX, 2L, object$data_info$X_scale,  "/")
    }
    X_new <- newX
    Z_new <- if (is.null(newZ)) X_new else {
      if (!is.matrix(newZ)) as.matrix(newZ) else newZ
    }
    on_train <- FALSE
  }

  if (is.null(t_eval)) t_eval <- object$breslow$times
  if (length(t_eval) == 0L)
    stop("No event times available in 'object$breslow$times'.")

  kernel_fn <- get_kernel(object$kernel)

  if (type %in% c("incidence", "survival", "all")) {
    u_new <- compute_index(X_new, object$gamma)
    p_new <- if (on_train) {
      object$p_hat
    } else {
      nw_predict(object$index_train, object$W,
                 u_new, object$h, kernel_fn)
    }
  }

  if (type %in% c("latency", "survival", "all")) {
    Su_new <- eval_Su(t_eval, object$breslow,
                      Z_new, object$beta)
  }

  if (type %in% c("survival", "all")) {
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

#' Compute the Prediction Error for the Incidence
#'
#' @description
#' Computes the out-of-sample prediction error for the
#' incidence part of a fitted SIC model on a test set,
#' as defined in Amico, Van Keilegom and Legrand (2019)
#' \doi{10.1111/biom.12999}:
#'
#' \deqn{PE = -\sum_{j=1}^{n_{\text{test}}} \left[
#'   \hat{w}_j \log \hat{p}(\mathbf{x}_j^{\text{test}}) +
#'   (1 - \hat{w}_j) \log\bigl(1 - \hat{p}(\mathbf{x}_j^{\text{test}})\bigr)
#' \right]}
#'
#' A smaller value indicates better prediction.
#'
#' @param object An object of class \code{"sic_fit"} fitted
#'   on the training set.
#' @param Y_test Numeric vector of follow-up times (test set).
#' @param delta_test Integer vector of event indicators (test set).
#' @param X_test Numeric matrix of incidence covariates (test set).
#' @param Z_test Numeric matrix of latency covariates (test set).
#'
#' @return A single numeric value: the prediction error (PE).
#'
#' @seealso \code{\link{predict.sic_fit}}, \code{\link{sic_fit}}
#'
#' @references
#' Amico M, Van Keilegom I, Legrand C (2019).
#' The single-index/Cox mixture cure model.
#' \emph{Biometrics}, \strong{75}(2), 452--462.
#' \doi{10.1111/biom.12999}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
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
#' fit <- sic_fit(Y[1:100], delta[1:100],
#'                X[1:100,], Z[1:100,,drop=FALSE],
#'                h_step = 0.1, verbose = 0)
#'
#' pe <- prediction_error(fit,
#'                        Y_test     = Y[101:150],
#'                        delta_test = delta[101:150],
#'                        X_test     = X[101:150,],
#'                        Z_test     = Z[101:150,,drop=FALSE])
#' cat("Prediction error:", round(pe, 3), "\n")
#' }
#'
#' @export
prediction_error <- function(object, Y_test, delta_test,
                             X_test, Z_test) {

  if (!is.matrix(X_test)) X_test <- as.matrix(X_test)
  if (!is.matrix(Z_test)) Z_test <- as.matrix(Z_test)

  p_test  <- predict.sic_fit(object, newX = X_test,
                              newZ = Z_test,
                              type = "incidence")
  Su_test <- eval_Su_diag(Y_test, Z_test,
                           object$beta, object$breslow)
  beyond  <- (delta_test == 0L) &
             (Y_test > object$data_info$tau)
  W_test  <- estep(p_test, Su_test, delta_test, beyond)
  p_test  <- clip_prob(p_test)
  -sum(W_test * safe_log(p_test) +
         (1.0 - W_test) * safe_log(1.0 - p_test))
}
