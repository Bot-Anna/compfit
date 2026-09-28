# Changelog

## compfit (development version)

### New features

- **[`reanchor_scenario()`](../reference/reanchor_scenario.md) extends a
  fitted model BACKWARDS in time.** It rewrites the fit’s sheet so the
  grid starts at an earlier `startpoint`, freezes every fitted parameter
  at its estimate (via [`fill_params()`](../reference/fill_params.md)),
  and pads the data with BLANK columns for the new years – blank cells
  are “missing” in the data-cell grammar, so nothing is invented and
  nothing enters the loss. The initial state now refers to the NEW start
  year, so respecify it with `states =` (a box is fitted, a bare number
  or expression is fixed); inheriting it warns, because the filled state
  describes the original start year. Returns a
  `modelParams`/`dataCombined` pair for
  [`fitCompartmentalModel()`](../reference/fitCompartmentalModel.md) or
  [`simulate_model()`](../reference/simulate_model.md).

  This is re-anchoring, not backward integration. The model is always
  solved FORWARD from the earlier year, so the earlier initial state is
  an assumption you supply or estimate. Integrating the ODE in reverse
  would be the literal reading of “project backwards”, but epidemic
  dynamics are contracting forwards and hence expanding in reverse, so a
  backward solve amplifies numerical error and leaves the feasible
  region (negative compartments, susceptibles above the population).
  Unlike `endpoint`, `startpoint` IS emitted into the generated ODE as a
  literal, which is why this returns a sheet to rebuild rather than
  reusing the fit.

- **The central curve of a Bayesian
  [`plot_fit()`](../reference/plot_fit.md) is now the posterior MEDIAN,
  and can be the pointwise median of the trajectory ensemble.**
  [`get_central_point()`](../reference/get_central_point.md) gains
  `summary = c("median", "mean")` and
  [`plot_fit()`](../reference/plot_fit.md) gains
  `central = c("median", "mean", "pointwise")`. Both previously used the
  mean only, which disagreed with `save_scenario()`’s own
  `summary = "median"` default and with the counterfactual reference;
  the median is also invariant under monotone reparameterisation and is
  not dragged into a tail by a skewed or funnel-shaped posterior (a
  `sigma`-to-zero funnel makes a mean meaningless).

  `"median"` and `"mean"` are PLUG-IN curves: solve the ODE once at that
  summary of the draws. A band is a different object – the 2.5%/97.5%
  quantiles of every draw’s output at each time point – and the two
  coincide only where the parameter-to-trajectory map is linear (means)
  or monotone (medians), which a compartmental model is not. So a
  plug-in line need not centre its own band. `central = "pointwise"`
  draws the ensemble’s pointwise median instead, which is that band’s
  50% quantile by construction. It costs an extra `n_draws` solves, and
  it is opt-in because a pointwise median is NOT a solution of the ODE:
  no single parameter vector produces it, so invariants such as a
  conserved `S+I+R` need not hold along it. Band frames now carry a
  `med` column alongside `lo`/`hi`.

  **Behaviour change:** an existing Bayesian
  [`plot_fit()`](../reference/plot_fit.md) call draws the median curve
  rather than the mean. For a well-identified posterior the difference
  is invisible; for a skewed one it is not. Pass `central = "mean"` for
  the old line.

- **Fits can be projected past their fitted horizon.**
  [`solve_and_evaluate()`](../reference/solve_and_evaluate.md) and
  [`plot_fit()`](../reference/plot_fit.md) take an `endpoint` argument:
  the model is solved to that year instead of the one on the sheet, so
  the trajectory – and, for a Bayesian fit, its bands or spaghetti lines
  – continues beyond the last observation while the data stays where it
  is. Extending the sheet instead is not an option, since a longer
  horizon there is rejected unless the data columns are padded to match.
  Nothing is recompiled: the generated ODE bakes in `startpoint` and
  `cutoff` but never `endpoint`, so the same registered model is solved
  over a longer grid, rebuilt exactly as the sheet would have built it.
  The years the two runs share are bit-identical. Truncating is rejected
  – the data would outrun the curve.

  This is projection, not forecasting: a time-varying `Function` keeps
  applying past the data (a ramp keeps ramping), so the extrapolation is
  only as meaningful as those formulas are outside the fitted window.

- **[`get_central_point()`](../reference/get_central_point.md) returns
  the central estimate of ANY successful fit** – the optimiser’s point
  estimate, or the posterior means for `method = "bayes"`. `fit$point`
  is `NULL` for a Bayesian fit, so the documented
  `solve_and_evaluate(fit, fit$point$initial_state, fit$point$parms)`
  idiom silently passed `NULL` and failed there. This is the vector
  [`plot_fit()`](../reference/plot_fit.md) already drew its central
  trajectory at; it is now reachable directly, and
  [`get_point()`](../reference/compfit-accessors.md)’s error message
  points to it.

- **The `Functions` column is sorted into dependency order.** Functions
  were emitted into the generated ODE verbatim, in sheet row order, and
  the block is straight-line code – so defining a Function below
  something that uses it failed. On the Julia backend that surfaced as
  `UndefVarError`; in the
  [`.derived_columns()`](../reference/dot-derived_columns.md) replay it
  was worse, because a failing entry is skipped by design, so the
  Function and everything downstream of it silently vanished from plots
  and data formulas. Row order now carries no meaning: the order is
  recovered from the reference graph instead, once, feeding the
  generated ODE and the replay alike, so all three backends agree. The
  sort is stable – a sheet already in a valid order generates
  byte-identical code – and permuting the column leaves the model
  numerically unchanged.

  Three faults cannot be fixed by reordering and are now rejected by
  [`validate_modelParams()`](../reference/validate_modelParams.md),
  before any backend runs: a name **defined more than once** (which
  definition a reference picks up would depend on row order), a
  **self-reference** (`x <- x + 1` – a Function is a definition, not an
  update), and a **reference cycle**, reported as the loop it forms
  (`'a' -> 'b' -> 'a'`) rather than as the full set of entries that
  could not be placed.

- **`dataCombined` accepts an `Asym` column.** `A+`/`A-` asymmetric data
  cells need a deviation, which previously had to be smuggled into the
  `Likelihood` column as `family; asym=<number>`. That overloaded a
  column named for the observation family with a penalty parameter, and
  it was read even under `method="mle"`, where the family token beside
  it is ignored. Put the number in its own `Asym` column instead: one
  positive value per stream. Precedence is a per-cell `A->B` deviation,
  then `Asym`, then the `asym=` clause – still honoured so existing
  sheets keep building, but no longer the documented form. The value is
  a **ratio of penalty slopes**, not a distance: `1` is symmetric, and
  `3` means a unit on the soft side costs a third of a unit on the hard
  side. Being a ratio, one value is right for a whole stream whatever
  its level. The metadata column names now live in a single
  `.DATA_META_COLS` constant, so a new meta column cannot be added in
  one place and forgotten in another (which would silently parse it as a
  year of observations).

### Breaking changes

- **Level columns are now `Level_<id>`, not `_`-prefixed.** A level
  (mixing subpopulation) column must start with `Level_`
  (e.g. `Level_1`, `Level_LA`, `Level_Age1`); the old bare-underscore
  form (`_Level1`, `_Age1`) is no longer detected. The id after `Level_`
  names the level and drives its companions (`N_<id>`, `Mixing_<id>`,
  `Pool_<id>`). Rename any `_<name>` level column to `Level_<name>`.

### Bug fixes

- **Mixing `[lo,hi]|init` and plain `[lo,hi]` fitted parameters no
  longer breaks post-fit recovery/plotting.** `.recover_solution()`
  required the fit’s fitted-parameter order to match `sap`’s exactly,
  but the Julia backend returns them in sheet order while `sap` groups
  with-init parameters first – so a sheet that interleaves the two
  tripped a “name/order mismatch” error in
  [`plot_fit()`](../reference/plot_fit.md) and friends. Fitted
  parameters are now reordered to `sap`’s canonical order by name (the
  p-vector the ODE indexes stays correct); only a genuine missing/extra
  fitted-parameter name is an error.

- **Censored / interval `negbin` (and `poisson`) data no longer break
  the Julia (NUTS) backend.** A censored or `[A,B]` cell on a
  discrete-family stream used `logcdf`/`logccdf`, which for
  `NegativeBinomial`/`Poisson` route through Rmath (a `Float64`-only C
  library with no ForwardDiff method) — so NUTS crashed with
  `MethodError: no method matching Float64(::ForwardDiff.Dual …)`. Those
  terms now use autodiff-safe pmf summation
  (`cf_disc_logcdf`/`cf_disc_logccdf`/ `cf_disc_loginterval`, built from
  `logpdf`, which is analytic and differentiable); intervals sum only
  over `(A,B]`. Continuous families (`gaussian`/`lognormal`) are
  unchanged (their CDF is already differentiable), and the MLE / R paths
  were never affected.

### New features

- **Initial-state expressions can reference state-independent
  `Functions`.** A `States` cell may now be written in terms of a
  `Functions` entry whose value is fixed by parameters and time alone –
  `sigma_*`, `q_*`, `sigmoid*`, `ntilde_*`, an equilibrium
  recent-fraction, etc. – e.g. `*R_H=(f_C+mu+sigma_HL)/(...)*Undiag`.
  The Function is evaluated at the initial time `t=-1` (calendar year
  `startpoint-1`), so time-varying ones collapse to their pre-window
  value (`q -> q_2014`, `sigmoid -> 0`). This removes the need to inline
  such expressions by hand. **State-DEPENDENT Functions** (`omega_*`,
  the `p_*` mixing probabilities, anything using a compartment or a
  level head `N_<level>`) remain unavailable at init – they are circular
  there – and
  [`validate_modelParams()`](../reference/validate_modelParams.md) now
  rejects a State cell that references one, naming the offending
  Function. The classification is automatic (a reference-graph fixpoint)
  and identical across the R, Julia, and Stan backends, which all build
  the same initial state.

- **Data/dummy formulas can reference the model’s `Functions` and level
  head counts.** A `dataCombined`/`dataDummy` formula is now evaluated
  on the solved trajectory in a scope that also contains every
  `Functions`-column quantity (e.g. `q_CH`, `omega_*`, the mixing
  probabilities `p_*`, the `sigmoid*`s) and the level head counts
  (`N1..Nk`, `N_<level>`, `total_pop`) – recomputed from the trajectory
  exactly as the ODE body computes them. So an overlay that needs a
  time-varying rate or a transmission denominator can name it
  (`annual(beta*H_world*q_CH*...)`) instead of inlining the definition.
  Quantities the replay cannot reconstruct are skipped individually, so
  the rest still resolve; a model with no `Functions` is unaffected.

- **State initial-value expressions may reference parameters by their
  plain name.** A parameter-dependent initial state can now be written
  `*S=Pi*Ntot` instead of requiring the `_0` alias form
  `*S=Pi_0*Ntot_0`; each declared name is rewritten to its `_0` alias
  internally, so both initial-state code paths resolve it. Existing
  `_0`-style sheets are unaffected.

- **Held-out initial-year data column.** `dataCombined` may now include
  a column named `startpoint - 1` (the initial-condition year, which is
  not a fitted snapshot). It is excluded from the fit but kept for
  plotting: [`plot_fit()`](../reference/plot_fit.md) draws its observed
  values as open markers, one year left of the first fitted point, so a
  pre-window data point can be shown for context without entering the
  likelihood. Absent such a column, nothing changes.

- **State names may be used directly in expressions.** Any
  `Linear`/`Quadratic`/ `Constant`/`Functions`/`Pool_`/`Mixing_` cell
  can now be written in terms of the compartments —
  e.g. `Hinf <- tau*R_HA + C_HA + delta*D_HA` in `Functions`, or a
  `Constant` inflow `0.02*R` — instead of positional `X[k]`. Each State
  name is rewritten to its slot in the generated model on every backend
  (R, Julia, Stan). A parameter/function may no longer share a State
  name ([`validate_modelParams()`](../reference/validate_modelParams.md)
  rejects it), since a compartment name always resolves to its state.

- **Level head counts are available by name.** Alongside the positional
  `N1..Nk`, each level now emits an alias `N_<level>` (e.g. `_Age1` →
  `N_Age1`), so an expression can refer to a stratum’s population by
  name instead of tracking its column order. `N_<level>` is reserved (a
  parameter/function may not take it).

- **`Mixing_<level>` columns** weight how each compartment enters a
  level’s transmission denominator. A companion to a `Level_<id>` column
  (e.g. `Mixing_1` for `Level_1`), it holds one weight per compartment
  defining an *effective* mixing pool `Nw = sum_k w_k X[k]` that
  normalises the second-order terms for that level, in place of the raw
  head count. A blank cell is the membership default (1 in-level, else
  0), so an absent or all-blank column reproduces the previous behaviour
  exactly. A non-zero weight on an out-of-level compartment pulls it
  into the pool (commuting / contact-matrix mixing). The raw `N<level>`
  and `total_pop` are left unchanged. Works on all three backends.

- **`Pool_<level>` columns** express a level’s transmission denominator
  as a whole *function of the level populations* – a single expression
  in `N1..Nk` / `total_pop` (plus parameters, `Functions`, and `time`),
  e.g. `N1 + c*N2/(1+N2/K)` for a saturating cross-level pool. It is the
  alternative to `Mixing_<level>` (a level uses one or the other) and is
  floored at a small positive value to guard against divide-by-zero;
  [`validate_modelParams()`](../reference/validate_modelParams.md) notes
  that the pool must stay positive across the solve. Level head counts
  `N1..Nk` and `total_pop` are now accepted symbols in coefficient /
  pool expressions.

### Breaking changes

- **Inline `$`-functions have been removed.** A coefficient cell (in
  `Linear<j>`, `Quadratic<j>`, or `Constant`) can no longer be a
  `$`-prefixed time-function (e.g. `$0.3*exp(-0.1*time)`). Define the
  expression as a named `Functions` entry and reference it by name
  instead – the general, reusable mechanism that already existed and
  works across every backend (the old `$` form only ever worked on the R
  backend).
  [`validate_modelParams()`](../reference/validate_modelParams.md) flags
  any remaining `$` cell with a migration hint. Consequently the codegen
  no longer generates the `f<ij>` / `g<ij>` / `cst<i>` helper names, and
  those are no longer reserved.

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
