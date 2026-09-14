# SICoxCURE: Single-Index/Cox Mixture Cure Model

<!-- badges: start -->
[![R-CMD-check](https://github.com/MohamedYassirErrahmani/SICoxCURE/actions/workflows/r.yml/badge.svg)](https://github.com/MohamedYassirErrahmani/SICoxCURE/actions/workflows/r.yml)
<!-- badges: end -->

## Overview

`SICoxCURE` implements the **Single-Index/Cox (SIC) mixture cure model**
proposed by Amico, Van Keilegom and Legrand (2019). The model writes
the population survival function as:

```
S(t | x, z) = 1 - p(x) + p(x) * Su(t | z)
```

where:
- `p(x) = g(γ'x)` is the **incidence** (probability of being uncured),
  modeled via a single-index structure with **unknown link function g**
  estimated nonparametrically (Nadaraya-Watson)
- `Su(t | z)` is the **latency** (survival of uncured subjects),
  modeled via a **Cox proportional hazards model**

Estimation uses an **EM algorithm** with bandwidth selected by
likelihood cross-validation at each iteration.

## Installation

```r
# From CRAN (once available)
install.packages("SICoxCURE")

# Development version from GitHub
# install.packages("devtools")
devtools::install_github("MohamedYassirErrahmani/SICoxCURE")
```

## Quick start

```r
library(SICoxCURE)

# Generate data
set.seed(42)
n     <- 200
X     <- matrix(rnorm(n * 4), n, 4)
Z     <- matrix(rbinom(n, 1L, 0.6), n, 1)
p     <- 1 / (1 + exp(-X %*% c(1, -0.5, 0.3, -0.2)))
B     <- rbinom(n, 1L, p)
T0    <- ifelse(B == 1, rexp(n, 0.5), Inf)
C     <- rexp(n, 0.3)
Y     <- pmin(T0, C)
delta <- as.integer(T0 <= C & T0 < Inf)

# Fit the SIC model
fit <- sic_fit(Y, delta, X, Z, verbose = 1)

# Summary with coefficients
summary(fit)

# Predictions on new data
pred <- predict(fit, newX = X[1:10, ], newZ = Z[1:10, ],
                type = "all")
```

## Functions

| Function | Description |
|----------|-------------|
| `sic_fit()` | Fit the SIC mixture cure model |
| `predict.sic_fit()` | Predict incidence, latency or survival |
| `prediction_error()` | Out-of-sample prediction error |
| `print.sic_fit()` | Compact display |
| `summary.sic_fit()` | Detailed summary with p-values |

## Reference

Amico M, Van Keilegom I, Legrand C (2019). The single-index/Cox
mixture cure model. *Biometrics*, **75**(2), 452–462.
https://doi.org/10.1111/biom.12999

## License

GPL-3 © 2025 Mohamed Yassir Errahmani, Mailis Amico,
Ingrid Van Keilegom, Catherine Legrand
