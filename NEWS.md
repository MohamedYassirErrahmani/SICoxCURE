# SICoxCURE 0.1.4

## New features

* **`constraint = "norm1"`** (new): implements the norm-1 identifiability constraint 
  (||gamma|| = 1) using `nloptr` with the SLSQP algorithm (equality constraint 
  applied at every EM iteration). Fallback to BFGS with re-normalization if nloptr 
  fails. More stable bootstrap standard errors than `gamma1` for all coefficients.

* **`constraint = "gamma1"`** (enhanced): now uses BFGS (as in Amico et al. 
  2019 reference code) instead of L-BFGS-B. Nelder-Mead fallback retained.

* **Brent method for bandwidth selection**: `h` is now selected by 
  `stats::optimize()` (Brent's method) at each EM iteration. Grid search is 
  retained as fallback if Brent fails. Tolerance controlled by `h_step`.

* **Formula interface**: `sic_fit(Surv(Y, d) ~ X1 + X2 | Z, data = dat)` is now 
  supported in addition to the matrix interface.

* **`max_iter_boot`**: new argument in `sic_fit()` and `bootstrap_se()` to control 
  the maximum number of EM iterations in each bootstrap replicate (default: `max_iter`).

* **`standardize_X = FALSE`** by default (previously `TRUE`). Set `standardize_X = TRUE` 
  to reproduce the standardization used in the original paper.

## Bug fixes and improvements

* License updated to 2026: Mohamed Yassir Errahmani & Maïlis Amico.
* `verbose` output improved: bandwidth grid, convergence criterion, and iteration 
  count are now displayed.
* `summary()` method improved: displays constraint type, standardization status, 
  and `call` information.

# SICoxCURE 0.1.3

## New features

* Bootstrap standard errors via `bootstrap_se()`.
* Parallel computing support via `parallel` argument.
* `plot_link()` and `plot_survival()` graphical functions.
* `predict.sic_fit()` for new observations.

# SICoxCURE 0.1.2

## Initial release

* `sic_fit()`: fit the SIC mixture cure model via EM algorithm.
* `estep()`: E-step of the EM algorithm.
* Nadaraya-Watson leave-one-out kernel estimator.
* Epanechnikov kernel (default) and Gaussian kernel.
