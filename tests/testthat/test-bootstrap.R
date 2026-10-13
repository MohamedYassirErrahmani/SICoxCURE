# =============================================================
# test-bootstrap.R
# Tests for bootstrap_se()
# =============================================================

make_data <- function(n = 80, seed = 99) {
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

test_that("bootstrap_se returns correct structure", {
  d   <- make_data(80, seed = 99)
  fit <- sic_fit(d$Y, d$delta, d$X, d$Z,
                 h_step = 0.1, verbose = 0)

  se  <- bootstrap_se(fit, B = 3L,
                      parallel = FALSE,
                      seed = 1L,
                      verbose = FALSE)

  expect_s3_class(se, "sic_bootstrap")
  expect_true(all(c("gamma","beta","gamma_boot",
                    "beta_boot","B","B_requested",
                    "n_failed") %in% names(se)))
  expect_length(se$gamma, length(fit$gamma))
  expect_length(se$beta,  length(fit$beta))
  expect_true(all(se$gamma >= 0))
  expect_true(all(se$beta  >= 0))
})

test_that("bootstrap_se SE are positive", {
  d   <- make_data(80, seed = 100)
  fit <- sic_fit(d$Y, d$delta, d$X, d$Z,
                 h_step = 0.1, verbose = 0)
  se  <- bootstrap_se(fit, B = 3L,
                      parallel = FALSE,
                      seed = 2L,
                      verbose = FALSE)
  expect_true(all(is.finite(se$gamma)))
  expect_true(all(is.finite(se$beta)))
})
