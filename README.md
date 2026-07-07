# compfit

<!-- badges: start -->
[![unit-tests](https://github.com/Bot-Anna/compfit/actions/workflows/tests.yml/badge.svg)](https://github.com/Bot-Anna/compfit/actions/workflows/tests.yml)
<!-- badges: end -->

**Fit, analyse and compare compartmental (SIR-type) ODE models from a
spreadsheet specification.**

`compfit` builds a compartmental ODE model from a small set of tables (initial
states, parameters, and per-compartment rate coefficients), fits it by **maximum
likelihood** or **Bayesian sampling**, and helps you analyse the result:
posterior summaries, prior-vs-posterior and predictive plots, global (Sobol) and
local sensitivity, practical identifiability, and self-contained code export.

Solves and Bayesian (NUTS) sampling run in **Julia** via the JuliaCall bridge for
speed — but Julia is **optional**: a pure-R backend (`deSolve`, with a
gradient-free sampler for the Bayesian case) fits without it, and a **Stan
backend** (`solver_control(backend = "stan")`, via `rstan`) offers gradient-based
NUTS with no Julia at all.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Bot-Anna/compfit")
```

The required R dependencies are installed automatically. To also pull the
optional packages used by particular backends and analyses, add
`dependencies = TRUE`. If the repository is private, set a token first with
`Sys.setenv(GITHUB_PAT = "<your-token>")` (a token with read access to the repo).

### Backends

- **Pure R** — no Julia, works out of the box: `solver_control(backend = "r")`.
- **Stan** — gradient-based NUTS, no Julia: `install.packages("rstan")`, then
  `solver_control(backend = "stan")`.
- **Julia** — fastest solves plus Julia/Turing NUTS. Install
  [Julia](https://julialang.org/) (>= 1.6), then once per session:

  ```r
  library(compfit)
  setup_julia()   # first call also installs the needed Julia packages
  ```

## Quick start

```r
library(compfit)

# A scenario is a folder with dataCombined / dataDummy / modelParams sheets:
sc <- load_scenario("path/to/MyScenario",
                    combined_file = "dataCombined.csv",
                    dummy_file    = "dataDummy.csv",
                    params_file   = "modelParams.csv")

# Maximum likelihood (works with or without Julia):
fit <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb")
summary(fit)
plot_fit(fit)$grid

# Bayesian (Julia/Turing NUTS):
fit_b <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "bayes",
                               bayes = bayes_control(chains = 4, iter = 2000))
posterior_report(fit_b)
plot_prior_posterior(fit_b)

# Bayesian without Julia: gradient-based NUTS via Stan (rstan):
fit_s <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "bayes",
                               solver = solver_control(backend = "stan"),
                               bayes  = bayes_control(chains = 4, iter = 2000))
```

### Bayesian backends

| backend | engine | needs | notes |
|---------|--------|-------|-------|
| `"julia"` (default) | Turing NUTS | Julia | fastest; full family/censoring/interval/asymmetric support |
| `"stan"` | Stan NUTS (`rstan`) | C++ toolchain | gradient-based, no Julia. All four families (`gaussian`/`lognormal`/`poisson`/`negbin`) and every data-cell type (observed, censored, interval `[A,B]`/`[A,B]~s`, asymmetric) — and it is the **strongest** backend for *censored counts* (Julia's autodiff can't handle those). |
| `"r"` | gradient-free MCMC (`BayesianTools`) | — | no Julia, no compiler; slower, mixes less well |

The same spreadsheet, priors, and post-fit tooling (`posterior_report`,
`plot_prior_posterior`, predictive plots) work across all three; only the sampler
differs. With `backend = "stan"` every R-side solve (plots, predictive) uses
deSolve, so no Julia session is ever needed.

Worked example scenarios (SI / SIS / SIR / SEIR and prior/family showcases) ship
with the package under `system.file("extdata", package = "compfit")`, each with a
README describing what it demonstrates.

## Simulate a fully-specified model (no data, no fitting)

If every state and parameter is fixed (`*name=value`), there is nothing to
estimate — you can run the model forward with **no `dataCombined` at all** and
evaluate any `dataDummy` formulas on the trajectory:

```r
sim_dir <- system.file("extdata", "SIR_sim", package = "compfit")
mp  <- read_data_file(file.path(sim_dir, "modelParams.csv"))
dd  <- read_data_file(file.path(sim_dir, "dataDummy.csv"))

sim <- simulate_model(mp, data_dummy = dd, solver = solver_control(backend = "r"))
sim                       # compartments, time span, fixed parameters
head(sim$sir_out)         # the solved trajectory
plot_simulation(sim)$grid # one panel per compartment / evaluated series
```

Fixed-twin `*_sim` scenarios (`minimal_sim`, `medium_sim`, `SI_sim`, `SIS_sim`,
`SIR_sim`, `SEIR_sim`) ship as ready-made examples. If any quantity is still
fittable, `simulate_model()` says so and points you at `fitCompartmentalModel()`.

## The model sheet, in one minute

`modelParams` uses column groups:

- `_Level*` — compartment column(s); reference each compartment by its 1-based
  index **or** by name (the two are interchangeable everywhere).
- `States` — initial values: `*S=val` (fixed), `S=[lo,hi]` (fitted), or
  `*S=<expr>` (parameter-dependent). **Compartments may be named freely** (`S`,
  `I`, `R`, ...) rather than `X1..Xn`; the **order of this column is canonical**
  (compartment 1 is the first entry). `X1..Xn` remains the identity case.
- `Parameters` — `*name=val` (fixed) or a **prior**: `name=[lo,hi]` (Uniform),
  `Normal(mu,sd)`, `LogNormal`, `Beta`, `Gamma`, `StudentT(nu,mu,sd)`, each
  optionally truncated with a trailing `[lo,hi]`.
- `Others` — `startpoint`, `endpoint`, `partition`, optional `cutoff`.
- `Functions` — time-varying helpers, e.g. `gamma_t<-gamma*(1+ramp*time)`.
- `Conditions` — constraints as a penalty, e.g. `beta>gamma`.
- `Linear<j>` / `Quadratic<j>` — first- and second-order rate coefficients, one
  pair of columns per compartment. `<j>` is the compartment index or its name
  (`LinearS` == `Linear1` when `S` is compartment 1). Blank cells are zero.

`dataCombined` has one row per observed stream (`Label`, `Formula`, and one
column per year), with an optional `Likelihood` family (`gaussian` / `poisson` /
`negbin` / `lognormal`) and censoring/interval/asymmetric data cells. The sheet
is checked at build time by `validate_modelParams()`, which reports entry
mistakes (non-numeric value, undeclared symbol in a coefficient, reversed box,
...) with a clear message naming the offending cell.

## Documentation

- `vignette("compfit")` — the full walkthrough.
- The PDF reference manual (`compfit_0.1.0.pdf`).
- `?fitCompartmentalModel`, `?plot_fit`, `?sobol_report`, `?identifiability_report`.

## License

GPL-3.

## Acknowledgments

`compfit` was developed with the assistance of [Claude](https://www.anthropic.com/claude),
Anthropic's AI assistant.
