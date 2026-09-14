# =============================================================
# test-predict.R
# Tests for predict.sic_fit() and prediction_error()
# =============================================================

make_data <- function(n = 100, seed = 10) {
  set.seed(seed)
  X     <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  Z     <- matrix(rbinom(n, 1L, 0.6), nrow = n, ncol = 1)
  p     <- 1 / (1 + exp(-X %*% c(1, -0.5, 0.3, -0.2)))
  B     <- rbinom(n, 1L, p)
  T_ev  <- rep(Inf, n)
  unc   <- which(B == 1L)
  if (length(unc) > 0)
    T_ev[unc] <- (-log(runif(length(unc))) /
                   (1.5 * exp(1.5 * Z[unc, 1])))^(1/1.2)
  C     <- rexp(n, 0.3)
  Y     <- pmin(T_ev, C)
  delta <- as.integer(T_ev <= C & T_ev < Inf)
  list(Y=Y, delta=delta, X=X, Z=Z)
}

test_that("predict returns correct types", {
  d   <- make_data(100, seed = 10)
  fit <- sic_fit(d$Y, d$delta, d$X, d$Z,
                 h_step = 0.1, verbose = 0)

  p <- predict(fit, type = "incidence")
  expect_true(is.numeric(p))
  expect_length(p, 100L)
  expect_true(all(p >= 0 & p <= 1))

  pred <- predict(fit, type = "all")
  expect_true(is.list(pred))
  expect_true(all(c("incidence","latency","survival",
                    "t_eval","index") %in% names(pred)))
})

test_that("predict on new data works", {
  d    <- make_data(100, seed = 11)
  fit  <- sic_fit(d$Y[1:80], d$delta[1:80],
                  d$X[1:80,], d$Z[1:80,,drop=FALSE],
                  h_step = 0.1, verbose = 0)
  pred <- predict(fit,
                  newX = d$X[81:100,],
                  newZ = d$Z[81:100,,drop=FALSE],
                  type = "incidence")
  expect_length(pred, 20L)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("prediction_error returns a positive scalar", {
  d   <- make_data(100, seed = 12)
  fit <- sic_fit(d$Y[1:80], d$delta[1:80],
                 d$X[1:80,], d$Z[1:80,,drop=FALSE],
                 h_step = 0.1, verbose = 0)
  pe  <- prediction_error(fit,
                          Y_test     = d$Y[81:100],
                          delta_test = d$delta[81:100],
                          X_test     = d$X[81:100,],
                          Z_test     = d$Z[81:100,,drop=FALSE])
  expect_true(is.numeric(pe))
  expect_length(pe, 1L)
  expect_true(pe > 0)
})
