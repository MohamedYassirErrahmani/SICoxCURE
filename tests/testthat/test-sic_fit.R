# =============================================================
# test-sic_fit.R
# Tests for sic_fit() -- using manual data generation
# =============================================================

# Helper: generate simple survival data manually
make_data <- function(n = 80, seed = 1) {
  set.seed(seed)
  X     <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  Z     <- matrix(rbinom(n, 1L, 0.6), nrow = n, ncol = 1)
  beta  <- 1.5
  gamma <- c(1, -0.5, 0.3, -0.2)
  p     <- 1 / (1 + exp(-X %*% gamma))
  B     <- rbinom(n, 1L, p)
  T_ev  <- rep(Inf, n)
  unc   <- which(B == 1L)
  if (length(unc) > 0)
    T_ev[unc] <- (-log(runif(length(unc))) /
                   (1.5 * exp(beta * Z[unc, 1])))^(1/1.2)
  C     <- rexp(n, 0.3)
  Y     <- pmin(T_ev, C)
  delta <- as.integer(T_ev <= C & T_ev < Inf)
  list(Y=Y, delta=delta, X=X, Z=Z)
}

test_that("sic_fit converges on simulated data", {
  d   <- make_data(80, seed = 1)
  fit <- sic_fit(d$Y, d$delta, d$X, d$Z,
                 h_step = 0.1, verbose = 0)
  expect_s3_class(fit, "sic_fit")
  expect_true(fit$converged)
  expect_length(fit$gamma, 4L)
  expect_length(fit$beta,  1L)
  expect_length(fit$p_hat, 80L)
  expect_true(all(fit$p_hat >= 0 & fit$p_hat <= 1))
})

test_that("sic_fit returns correct structure", {
  d   <- make_data(80, seed = 2)
  fit <- sic_fit(d$Y, d$delta, d$X, d$Z,
                 h_step = 0.1, verbose = 0)
  expect_true(all(c("gamma","beta","g_hat","p_hat","Su_hat",
                    "breslow","h","W","loglik","converged",
                    "n_iter","crit") %in% names(fit)))
  expect_true(fit$h >= 0.4 && fit$h <= 1.0)
  expect_true(all(fit$W >= 0 & fit$W <= 1))
  expect_true(all(fit$W[d$delta == 1] == 1.0))
})

test_that("sic_fit data.frame interface works", {
  d   <- make_data(80, seed = 3)
  dat <- data.frame(Y=d$Y, delta=d$delta,
                    X1=d$X[,1], X2=d$X[,2],
                    X3=d$X[,3], X4=d$X[,4],
                    Z=d$Z[,1])
  fit <- sic_fit(data      = dat,
                 col_Y     = "Y",
                 col_delta = "delta",
                 cols_X    = c("X1","X2","X3","X4"),
                 cols_Z    = "Z",
                 h_step    = 0.1,
                 verbose   = 0)
  expect_s3_class(fit, "sic_fit")
  expect_true(fit$converged)
})

test_that("sic_fit with constraint norm1 works", {
  d   <- make_data(80, seed = 4)
  fit <- sic_fit(d$Y, d$delta, d$X, d$Z,
                 constraint = "norm1",
                 h_step = 0.1, verbose = 0)
  expect_equal(sqrt(sum(fit$gamma^2)), 1.0, tolerance = 1e-6)
  expect_true(fit$gamma[1] > 0)
})

test_that("sic_fit missing arguments raise errors", {
  expect_error(sic_fit(), "missing")
  expect_error(sic_fit(Y = c(1,2,3)), "missing")
})
