# Changelog

## compfit 0.1.0

First release.

### Features

- Build compartmental (SIR-type) ODE models from a spreadsheet
  specification and fit them by maximum likelihood (`lbfgsb` / `deoptim`
  / `hypercube`) or Bayesian sampling
  ([`fitCompartmentalModel()`](../reference/fitCompartmentalModel.md)).
- Analysis tools: posterior summaries and diagnostics
  ([`posterior_report()`](../reference/posterior_report.md),
  [`posterior_draws()`](../reference/posterior_draws.md), …),
  prior-vs-posterior and predictive plots
  ([`plot_prior_posterior()`](../reference/plot_prior_posterior.md),
  [`plot_fit()`](../reference/plot_fit.md)), global (Sobol) and local
  sensitivity ([`sobol_report()`](../reference/sobol_report.md),
  [`local_sensitivity()`](../reference/local_sensitivity.md)), practical
  identifiability
  ([`identifiability_report()`](../reference/identifiability_report.md)),
  counterfactual scenarios (`run_counterfactual_folder()` and friends),
  and self-contained code export
  ([`extract_code()`](../reference/extract_code.md)).

### Solver backends

- `solver_control(backend = ...)` selects the ODE/inference backend:
  - `"julia"` (default) — fast solves via OrdinaryDiffEq and
    gradient-based (NUTS) Bayesian sampling via Turing, through the
    JuliaCall bridge. Requires a Julia runtime; initialise with
    [`setup_julia()`](../reference/setup_julia.md).
  - `"r"` — a pure-R backend (deSolve) that needs **no Julia**. Maximum
    likelihood is fully supported. Bayesian sampling uses a
    gradient-free MCMC (BayesianTools DEzs), seeded from a quick MLE
    fit; it is **much slower and mixes far worse than the Julia/NUTS
    path** and is intended for small models or quick checks only.
- `JuliaCall` and `diffeqr` are now `Suggests`: the package installs and
  runs the R backend without them.
  [`setup_julia()`](../reference/setup_julia.md) reports clearly if they
  (or Julia) are missing.

### Notes

- Example scenarios ship as package data under `inst/extdata/`
  (reachable via
  `system.file("extdata", "minimal", package = "compfit")`).
- [`extract_code()`](../reference/extract_code.md) emits human-readable
  data by default (`inline = "readable"`); `inline = "fidelity"` keeps
  the byte-exact [`serialize()`](https://rdrr.io/r/base/serialize.html)
  form.
