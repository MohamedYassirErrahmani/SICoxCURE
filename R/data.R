#' Simulated SIC Mixture Cure Dataset
#'
#' A simulated dataset illustrating the use of the SICoxCURE package.
#' Data were generated from a Single-Index/Cox mixture cure model with
#' n = 300 observations, three incidence covariates (TRT, SEX, AGE),
#' and one latency covariate (TRT).
#'
#' @format A data frame with 300 rows and 5 variables:
#' \describe{
#'   \item{FAILTIME}{Observed follow-up time (numeric)}
#'   \item{FAILCENS}{Event indicator: 1 = event, 0 = censored (integer)}
#'   \item{TRT}{Treatment indicator: 0 or 1 (integer)}
#'   \item{SEX}{Sex indicator: 0 or 1 (integer)}
#'   \item{AGE}{Age in years (numeric)}
#' }
#' @details
#' Simulated from a SIC mixture cure model with:
#' \itemize{
#'   \item Incidence: logit-1(gamma'x - 0.5) with gamma = (0.8, -0.3, 0.015).
#'   \item Latency: exponential with rate 0.25 * exp(-0.6 * TRT).
#'   \item Censoring: exponential with rate 0.12.
#'   \item Event rate: 23.3%; plateau: 18.7%; seed: 2019.
#' }
#' @examples
#' \donttest{
#' data(sicure_data)
#' fit <- sic_fit(
#'   Y             = sicure_data$FAILTIME,
#'   delta         = sicure_data$FAILCENS,
#'   X             = as.matrix(sicure_data[, c('TRT','SEX','AGE')]),
#'   Z             = matrix(sicure_data$TRT, ncol = 1),
#'   standardize_X = TRUE,
#'   verbose       = 1L
#' )
#' print(fit)
#' }
"sicure_data"
