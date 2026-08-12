# compfit (development version)

## Breaking changes

* **Level columns are now `Level_<id>`, not `_`-prefixed.** A level (mixing
  subpopulation) column must start with `Level_` (e.g. `Level_1`, `Level_LA`,
  `Level_Age1`); the old bare-underscore form (`_Level1`, `_Age1`) is no longer
  detected. The id after `Level_` names the level and drives its companions
  (`N_<id>`, `Mixing_<id>`, `Pool_<id>`). Rename any `_<name>` level column to
  `Level_<name>`.

## New features

* **Held-out initial-year data column.** `dataCombined` may now include a column
  named `startpoint - 1` (the initial-condition year, which is not a fitted
  snapshot). It is excluded from the fit but kept for plotting: `plot_fit()` draws
  its observed values as open markers, one year left of the first fitted point, so
  a pre-window data point can be shown for context without entering the likelihood.
  Absent such a column, nothing changes.
* **State names may be used directly in expressions.** Any `Linear`/`Quadratic`/
  `Constant`/`Functions`/`Pool_`/`Mixing_` cell can now be written in terms of the
  compartments — e.g. `Hinf <- tau*R_HA + C_HA + delta*D_HA` in `Functions`, or a
  `Constant` inflow `0.02*R` — instead of positional `X[k]`. Each State name is
  rewritten to its slot in the generated model on every backend (R, Julia, Stan).
  A parameter/function may no longer share a State name (`validate_modelParams()`
  rejects it), since a compartment name always resolves to its state.
* **Level head counts are available by name.** Alongside the positional `N1..Nk`,
  each level now emits an alias `N_<level>` (e.g. `_Age1` → `N_Age1`), so an
  expression can refer to a stratum's population by name instead of tracking its
  column order. `N_<level>` is reserved (a parameter/function may not take it).

* **`Mixing_<level>` columns** weight how each compartment enters a level's
  transmission denominator. A companion to a `Level_<id>` column (e.g.
  `Mixing_1` for `Level_1`), it holds one weight per compartment defining an
  *effective* mixing pool `Nw = sum_k w_k X[k]` that normalises the second-order
  terms for that level, in place of the raw head count. A blank cell is the
  membership default (1 in-level, else 0), so an absent or all-blank column
  reproduces the previous behaviour exactly. A non-zero weight on an out-of-level
  compartment pulls it into the pool (commuting / contact-matrix mixing). The raw
  `N<level>` and `total_pop` are left unchanged. Works on all three backends.
* **`Pool_<level>` columns** express a level's transmission denominator as a whole
  *function of the level populations* -- a single expression in `N1..Nk` /
  `total_pop` (plus parameters, `Functions`, and `time`), e.g. `N1 + c*N2/(1+N2/K)`
  for a saturating cross-level pool. It is the alternative to `Mixing_<level>` (a
  level uses one or the other) and is floored at a small positive value to guard
  against divide-by-zero; `validate_modelParams()` notes that the pool must stay
  positive across the solve. Level head counts `N1..Nk` and `total_pop` are now
  accepted symbols in coefficient / pool expressions.

## Breaking changes

* **Inline `$`-functions have been removed.** A coefficient cell (in `Linear<j>`,
  `Quadratic<j>`, or `Constant`) can no longer be a `$`-prefixed time-function
  (e.g. `$0.3*exp(-0.1*time)`). Define the expression as a named `Functions` entry
  and reference it by name instead -- the general, reusable mechanism that already
  existed and works across every backend (the old `$` form only ever worked on the
  R backend). `validate_modelParams()` flags any remaining `$` cell with a
  migration hint. Consequently the codegen no longer generates the `f<ij>` /
  `g<ij>` / `cst<i>` helper names, and those are no longer reserved.

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
