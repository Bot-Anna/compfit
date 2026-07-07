------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# Test suite — companion index

The `compfit` package test suite, on **testthat** (edition 3). A small helper
(`tests/testthat/helper-compfit.R`) re-implements the original harness API
(`chk`, `chk_equal`, `chk_error`, `chk_ok`, the `th_*` flow helpers) on top of
testthat expectations, so each test file reads the same as before. Pure-R unit
tests run anywhere; the integration tests start Julia via `setup_julia()` and
**skip cleanly** when Julia (or a scenario) is unavailable.

## How to run

``` r
# Everything (integration tests start Julia; skip if unavailable):
devtools::test()

# Just the pure-R unit tests (fast, no Julia):
devtools::test(filter = "unit")

# A single file:
devtools::test(filter = "unit-utils")

# Or the full check:
devtools::check()
```

``` bash
# Point the integration tests at your own scenario instead of the fixture:
FITCM_SCENARIO_DIR="/path/to/Scenario" Rscript -e 'devtools::test()'
```

The **interactive workflow is unchanged**: `source("setup.R")` (after setting
`PROJECT_ROOT`) still loads everything and initialises Julia.

## What each test file checks

| File | Needs | What it verifies |
|----|----|----|
| `testthat/helper-compfit.R` | — | Not a test: the testthat-backed harness shim (`chk`, `chk_equal`, `chk_error`, `chk_ok`, `th_skip`→`skip()`), project-root discovery, and `th_load_pure()` (no-op; package preloaded) / `th_have_setup()` (`setup_julia()`). |
| `test-unit-utils.R` | base R | `parse_data_cell` grammar (observed / missing / `<L` `<=L` `>L` `>=L` censoring, whitespace, unparseable error); `parseLikelihood` families + aliases + dispersion flag; `parsePrior` (box / Normal / truncated / fixed); family-registry constants; `normalise`/`denormalise`, `is_cumulative_stream`, `extract_param`, `extract_numbers`, `snap_to_step`, `remove_trailing_plus`, `reduce_expression`. |
| `test-unit-scenario.R` | base R | `next_plot_path` numbering (first→`_1`, max+1, legacy un-numbered ignored, regex-special stems); `read_data_file(text_cols=)` preserves stored-as-text numbers (the intermittent-`0` fix) and round-trips via `as.numeric`; unsupported-extension error; `save_grid` no-op when `$grid` is `NULL`. |
| `test-unit-modelbuild.R` | base R | `numberOfComps` coerces `_Level` indices numerically (so `"10" > "2"`, not lexicographic); `.bounds` for named fitted quantities (midpoint `init_norm`), the **fully-fixed sheet** case (zero fitted → empty bounds, no error), and the unnamed-bounds error. |
| `test-unit-prepare-data.R` | base R | `.prepare_data`: left/right censoring matrices (`cens_mask`/`lcens_mask`/`limit_mat`/`llimit_mat`/`inc_mask`/`linc_mask`), years×streams orientation, cumulative differencing + `cumulative_cols`, `Weight`/`Average` numeric coercion (text `"0"`→0; non-numeric→0 **with warning**), and the column-count-vs-horizon guard. |
| `test-unit-evaluate.R` | base R | `evaluate_formula` on a synthetic trajectory: bare stocks, bare parameter names, arithmetic, `annual()` rolling integral (incl. leading `NA`s), `cumulative()` running integral. |
| `test-unit-recover-solution.R` | base R | `.recover_solution` rebuilds natural-scale `(initial_state, parms)`; the order-contract guard trips on a name/order mismatch but is empty-safe for a fully-fixed (counterfactual) sheet where all quantities are fixed. || `test-unit-fillparams.R` | base R (+ `writexl`,`readxl` for the last block) | `.fp_rewrite_cell` (fitted→`*name=value`, fixed/unknown/NA untouched, digit rounding); `fill_params` rewrites States/Parameters and leaves structural columns intact; **optional** `write_filled_params` + `verify_filled` xlsx round-trip. |
| `test-unit-extract-code.R` | base R | Guards `extract_code`'s inlined-objects list (`.xc_loss_captured`) against `lossFunction`'s formals — fails if a captured data argument is added to the loss but not to the emitted script (the class of bug that hid `lcens_mask`/`llimit_mat`). |
| `test-unit-local-sensitivity.R` | base R | `.ls_split_point` overlays a perturbed natural-scale vector onto a base point (replaces fitted params/states, keeps fixed ones; no-op when nothing overlaps). |
| `test-unit-sobol-report.R` | base R | `sobol_report` reporting logic on a synthetic tidy Sobol table — validity (in-`[0,1]`), the interaction column, global parameter ranking, not-converged flagging for out-of-range indices, and that a precomputed data.frame is reported without recomputing. |
| `test-unit-modelchain.R` | base R (+ `readr` optional) | Builds the minimal fixture through the whole parsing/codegen chain (`numberOfComps` → `statesAndParams` → `generateExpressions` → `compartmentalFunction`) **without** Julia, asserting the emitted Julia ODE source has both derivatives and the linear coefficients. Validates the fixture and a large slice of model-building code. |
| `test-unit-modelchain-medium.R` | base R (+ `readr` optional) | Builds the **medium** fixture (3-compartment SIR) through the chain, exercising what the minimal model lacks: a quadratic infection term, a time-varying recovery `Function`, a `Condition`→penalty, and **parameter-dependent function-defined initial states** (asserts the state functions compute `X1=N0·(1−init_inf)`, `X2=init_inf·N0`). |
| `test-unit-datacell-grammar.R` | base R (+ `deSolve`) | Interval `[A,B]`, **soft** interval `[A,B]~s`, and asymmetric `A->B`/`A+`/`A-` data cells: `.prepare_data` masks/edges/shoulder scales, an MLE fit on the R backend with mixed cells, the `A+`/`A-` "no `asym=` declared" error, no-op equivalences with right-censoring, `plot_fit` rendering, and the generated Turing `@model` string carrying the new branches (incl. the no-CDF soft-interval branch). |
| `test-unit-coverage.R` | base R (+ `deSolve`) | The coverage/scaling pattern (`observed = rho * true` via a `Formula` like `rho*X2`): `rho` is routed to the fitted parameters (not the param-functions), and an MLE fit with an anchoring stream recovers `rho`. |
| `test-unit-validate.R` | base R | `validate_modelParams()`: every shipped fixture passes (no false positives), and each class of bad entry raises a clear, specific error — non-numeric fixed value / box / prior arg / time-grid field, reversed box, out-of-range `Quadratic` target, undefined symbol in a coefficient, missing time-grid key. |
| `test-unit-examples.R` | base R | Guards the four shipped epidemic example scenarios (`SI`/`SIS`/`SIR`/`SEIR` under `inst/extdata/`): each `load_scenario` → `.prepare_data` succeeds, the compartment count is right, and the special data cells each example advertises (`<=`/`>=`/`[A,B]`/`A->B`/`x`) and its declared families land in the expected masks. |
| `test-integration-fit.R` | **Julia + scenario** | `load_scenario` → MLE `fitCompartmentalModel` (class/success/point/censoring masks) → `save_fit`/`load_fit` round-trip (loss rebuilt, error preserved) → `plot_fit` (plots/code/size) → `extract_code`. |
| `test-integration-datacell-grammar.R` | **Julia** | Real Turing/Julia round-trip for the interval/asymmetric grammar: a Bayesian fit with mixed observed/interval/asym cells registers the generated `@model` and samples it; `model_code` carries the `cf_logsubexp`/`asym_dir_mat` branches. |
| `test-integration-medium.R` | **Julia** | Fits the medium SIR fixture end to end; checks the function-defined initial states resolve (`X1=990`, `X2=10`, `X3=0`), the trajectory is finite/non-negative, and **S+I+R is conserved** (the quadratic flux routing is physically correct). |
| `test-integration-sensitivity.R` | **Julia** (+ `sensitivity` for the Sobol part) | `local_sensitivity` at the fit (structure, tornado/trajectory plots) on the medium SIR, including the physical check that the conserved total `S+I+R` has ~zero sensitivity; plus a small `sobol_loss` run returning first-order/total indices. |
| `test-integration-extract-code.R` | **Julia + scenario** | Emits the self-contained `extract_code(what="all")` script, runs it in a fresh environment, and confirms it re-registers Julia, re-fits, and **reproduces the fitted quantities** within tolerance. |
## Configuration (environment variables)

| Variable | Used by | Meaning |
|----|----|----|
| `FITCM_ROOT` | all | Project root, if not running from it. |
| `FITCM_SCENARIO_DIR` | integration | Folder with the input workbooks. **Defaults to the committed minimal fixture** (`inst/extdata/minimal/`, resolved via `system.file()`), so the integration tests run end-to-end wherever Julia works. Set this to point at your own scenario instead. |
| `FITCM_COMBINED` / `FITCM_DUMMY` / `FITCM_PARAMS` | integration | Override input file names (defaults: `dataCombined.xlsx`, `dataDummy.xlsx`, `model.xlsx` → `modelParams.xlsx`). |
## Expected workbooks (for the integration tests)

In `FITCM_SCENARIO_DIR`: - `dataCombined.xlsx` (or your `FITCM_COMBINED`) — required - model parameter sheet `model.xlsx` / `modelParams.xlsx` (or `FITCM_PARAMS`) — required - `dataDummy.xlsx` (or `FITCM_DUMMY`) — optional; some checks skip without it

Two committed fixtures ship as package data under `inst/extdata/` (CSV, so
they're git-friendly), reachable from installed code via
`system.file("extdata", "minimal", package = "compfit")`:

- `minimal/` — a 2-compartment **linear** model. The default for `FITCM_SCENARIO_DIR`; covers the basic pipeline.
- `medium/` — a 3-compartment **SIR** building on the minimal one: a quadratic infection term, a time-varying recovery `Function`, a `Condition` penalty, and parameter-dependent function-defined initial states. Used by `*-medium.R`.

Four worked **epidemic example scenarios** (fit by both MLE and Bayes; each has
its own `README.md`), guarded by `test-unit-examples.R`:

- `SI/` — 2-compartment mass-action infection; **poisson** counts, a missing cell.
- `SIS/` — S↔I with a time-varying recovery `Function` and a `Condition`; **gaussian**, with `<=` and `A->B` cells.
- `SIR/` — two streams / two families (**negbin** `I` + **gaussian** `R`); interval `[A,B]` and `>=` cells on the continuous `R` stream.
- `SEIR/` — 4-compartment chain, three families (**gaussian**/**negbin**/**lognormal**); `>=`, `[A,B]`, and missing cells.
- `SIR_priors/` — SIR + reporting coverage; prior showcase (`StudentT`/`Normal`/`Beta`/`Uniform`); **negbin** + **gaussian**.
- `SEIR_priors/` — SEIR; prior showcase (`Uniform`/`Gamma`/`LogNormal`/`Beta`); **gaussian** + **poisson** + **lognormal**. Together with `SIR_priors` covers every prior distribution and every likelihood family.

  These deliberately keep CDF-based cells (`<=`/`>=`/`[A,B]`/`A->B`) on **continuous**-family streams: the Julia/NUTS path cannot autodiff `logcdf` for discrete families, so `poisson`/`negbin` streams carry plain observed/missing cells only (the R backend has no such restriction).

Both run without any private data — only a working Julia bridge is required for the integration parts.

## Notes

- Unit tests are pure R and run anywhere — start with `devtools::test(filter = "unit")`.
- A skip is not a failure; it means a prerequisite (Julia, a workbook, an optional package) was absent.
- Integration tests perform real ODE solves and can take a while (especially any Bayesian step).
- Continuous integration: `.github/workflows/tests.yml` runs the pure-R unit tests on every push / PR (no Julia needed); the integration tests auto-skip there.
