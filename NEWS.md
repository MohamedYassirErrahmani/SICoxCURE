# SICoxCURE 0.1.0

## First release

* Implements the Single-Index/Cox (SIC) mixture cure model of
  Amico, Van Keilegom and Legrand (2019).

* `sic_fit()`: main function to fit the SIC model via EM algorithm.
  Accepts both vector/matrix and data.frame interfaces.

* `predict.sic_fit()`: predictions for incidence, latency and
  population survival on new observations.

* `prediction_error()`: out-of-sample prediction error for the
  incidence (as in Article Section 4).

* `get_link_function()`: true link function for each scenario.

* `compute_ase()`: Average Squared Error for evaluating incidence
  estimation accuracy.

## Improvements over the original code (Amico et al. 2019)

* `init_latency()`: corrected to fit Cox model on uncensored
  observations only, as specified in the article (p.456).

* `select_bandwidth()`: replaced fixed grid of 30 values with a
  user-controlled step size (`h_step = 0.001` by default, giving
  601 candidate values in [0.4, 1.0]).

* `mstep_incidence()` and `mstep_latency()`: added L-BFGS-B
  optimisation with Nelder-Mead fallback for robustness.

* `convergence_criterion()`: added `method` argument
  ("both", "max", "sum") for flexibility.
