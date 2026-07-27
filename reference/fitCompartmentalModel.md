# Fit a compartmental model

Build a compartmental ODE model from a spreadsheet specification and fit
it by maximum likelihood or Bayesian sampling. By default ODE solves
(and Julia NUTS) run in Julia; call \[setup_julia()\] first or rely on
the lazy initialiser. A pure-R backend and a Stan backend need no Julia
– see \[solver_control()\].

## Usage

``` r
fitCompartmentalModel(
  modelParams,
  dataCombined,
  method = c("lbfgsb", "deoptim", "hypercube", "bayes"),
  solver = solver_control(),
  control = optim_control(),
  hypercube = hypercube_control(),
  bayes = NULL,
  init = NULL,
  checkpoint_file = "checkpoint_best_solution.rds",
  debug_env = NULL
)
```

## Arguments

- modelParams:

  Model-parameter sheet (data frame): compartments, priors, rate
  coefficients and the time grid. See \`vignette("compfit")\` for the
  grammar.

- dataCombined:

  Observed data, one row per stream (\`Formula\` + one column per year,
  optional \`Label\`/\`Likelihood\`/\`Weight\`/\`Average\`). Required
  for fitting; to run a fully-fixed model forward with no data use
  \[simulate_model()\] instead.

- method:

  Fitting method (one of):

  \`"lbfgsb"\`

  :   deterministic local MLE (L-BFGS-B). Fast; the default.

  \`"deoptim"\`

  :   stochastic global MLE (differential evolution). More robust on
      rugged/multi-modal objectives, slower. Seed via
      \[optim_control()\].

  \`"hypercube"\`

  :   deterministic Sobol pre-search; a coarse global scan (see
      \[hypercube_control()\]).

  \`"bayes"\`

  :   full Bayesian posterior sampling (see \[bayes_control()\] and the
      \`backend\` in \[solver_control()\]).

- solver:

  Solver/backend settings, see \[solver_control()\] – chooses
  \`"julia"\` / \`"r"\` / \`"stan"\` and the ODE tolerances.

- control:

  Optimiser settings for the MLE methods, see \[optim_control()\].

- hypercube:

  Pre-search settings for \`method = "hypercube"\`, see
  \[hypercube_control()\].

- bayes:

  Sampler settings, see \[bayes_control()\]; used only for \`method =
  "bayes"\` (defaults are supplied if \`NULL\`).

- init:

  Optional starting point as a named numeric vector on the normalised
  \\0,1\\ scale (advanced; normally left \`NULL\`).

- checkpoint_file:

  Path to an \`.rds\` where the running best solution is written during
  an MLE search, so a long fit can be recovered if interrupted.

- debug_env:

  Optional environment into which intermediate build objects are stashed
  for debugging (no effect on the result when \`NULL\`).

## Value

An object of class \`"compartmentalFit"\`: \`\$point\` (natural-scale
estimates) for MLE, \`\$samples\` for Bayes, plus the model, data,
bounds and solver used. Inspect with \[summary()\], \[get_point()\],
\[posterior_report()\], \[plot_fit()\].

## Examples

``` r
mini <- system.file("extdata", "minimal", package = "compfit")
sc   <- load_scenario(mini, combined_file = "dataCombined.csv",
                      dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
# Julia-free maximum-likelihood fit (pure-R deSolve backend) -- runs at check:
fit <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb",
                             solver = solver_control(backend = "r"),
                             checkpoint_file = tempfile(fileext = ".rds"))
#> [1] 16555
#> [1] 16538.7
#> [1] 16523.32
#> [1] 4477.895
#> [1] 4466.666
#> [1] 4451.289
#> [1] 1424.215
#> [1] 1421.278
#> [1] 1418.881
#> [1] 1411.398
#> [1] 1342.575
#> [1] 1338.998
#> [1] 1338.657
#> [1] 1331.189
#> [1] 1089.118
#> [1] 1086.611
#> [1] 1085.415
#> [1] 1079.285
#> [1] 758.9525
#> [1] 752.2405
#> [1] 747.0081
#> [1] 746.6348
#> [1] 545.2051
#> [1] 543.9542
#> [1] 540.9125
#> [1] 422.513
#> [1] 419.3855
#> [1] 365.6516
#> [1] 362.9373
#> [1] 362.8489
#> [1] 264.1457
#> [1] 263.3836
#> [1] 262.5371
#> [1] 258.3689
#> [1] 258.1481
#> [1] 257.3876
#> [1] 256.4764
#> [1] 256.4059
#> [1] 256.3689
#> [1] 256.355
#> [1] 256.3548
#> [1] 256.3548
get_point(fit)
#> $initial_state
#>       X1       X2 
#> 115.0665   0.0000 
#> 
#> $parms
#>         k         m 
#> 0.1782059 0.0191363 
#> 
if (FALSE) { # \dontrun{
# The Julia backend: faster, and required for NUTS Bayesian sampling.
setup_julia()
fit_b <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "bayes",
                               bayes = bayes_control(chains = 2, iter = 1000))
} # }
```
