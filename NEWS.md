# compfit 0.1.0

First release.

## Features

* Build compartmental (SIR-type) ODE models from a spreadsheet specification and
  fit them by maximum likelihood (`lbfgsb` / `deoptim` / `hypercube`) or Bayesian
  sampling (`fitCompartmentalModel()`).
* Analysis tools: posterior summaries and diagnostics (`posterior_report()`,
  `posterior_draws()`, ...), prior-vs-posterior and predictive plots
  (`plot_prior_posterior()`, `plot_fit()`), global (Sobol) and local sensitivity
  (`sobol_report()`, `local_sensitivity()`), practical identifiability
  (`identifiability_report()`), counterfactual scenarios
  (`run_counterfactual_folder()` and friends), and self-contained code export
  (`extract_code()`).

## Solver backends

* `solver_control(backend = ...)` selects the ODE/inference backend:
  * `"julia"` (default) --- fast solves via OrdinaryDiffEq and gradient-based
    (NUTS) Bayesian sampling via Turing, through the JuliaCall bridge. Requires a
    Julia runtime; initialise with `setup_julia()`.
  * `"r"` --- a pure-R backend (deSolve) that needs **no Julia**. Maximum
    likelihood is fully supported. Bayesian sampling uses a gradient-free MCMC
    (BayesianTools DEzs), seeded from a quick MLE fit; it is **much slower and
    mixes far worse than the Julia/NUTS path** and is intended for small models
    or quick checks only.
* `JuliaCall` and `diffeqr` are now `Suggests`: the package installs and runs the
  R backend without them. `setup_julia()` reports clearly if they (or Julia) are
  missing.

## Notes

* Example scenarios ship as package data under `inst/extdata/` (reachable via
  `system.file("extdata", "minimal", package = "compfit")`).
* `extract_code()` emits human-readable data by default (`inline = "readable"`);
  `inline = "fidelity"` keeps the byte-exact `serialize()` form.
