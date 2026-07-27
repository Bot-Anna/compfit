# Bayesian sampling settings

Bundle the settings for \`method = "bayes"\`. These apply across all
three Bayesian backends selected by \[solver_control()\]: Julia/Turing
NUTS (\`backend = "julia"\`), Stan/rstan NUTS (\`backend = "stan"\`),
and the gradient-free \`BayesianTools\` sampler (\`backend = "r"\`).
Priors themselves are read from the model sheet's
\`Parameters\`/\`States\` columns (per-quantity boxes or named
distributions); the settings here govern the sampler, not the priors.

## Usage

``` r
bayes_control(
  sampler = "NUTS(0.65)",
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = NULL,
  init_from_optim = TRUE,
  progress = TRUE,
  sigma_prior = "truncated(Normal(0, 1), 0, Inf)",
  phi_prior = "Gamma(2, 5)",
  sigma_prior_stan = "normal(0, 1)",
  phi_prior_stan = "gamma(2, 0.2)",
  adapt_delta = 0.8
)
```

## Arguments

- sampler:

  Julia/Turing sampler expression, e.g. \`"NUTS(0.65)"\` (the argument
  is the target acceptance rate). Used by \`backend = "julia"\` only;
  Stan uses its own NUTS and the R backend uses DEzs.

- chains:

  Number of independent MCMC chains (default 4). Multiple chains enable
  the R-hat convergence diagnostic.

- iter:

  Total iterations per chain, \*\*including\*\* warmup (default 2000).

- warmup:

  Warmup / burn-in iterations per chain, discarded before inference
  (default 1000); the kept sample size per chain is \`iter - warmup\`.

- seed:

  Optional RNG seed for reproducible sampling.

- init_from_optim:

  If \`TRUE\` (default), run a quick MLE fit first and initialise the
  sampler there (a MAP-style warm start). Applies to the Julia and R
  backends; the Stan backend relies on Stan's own warmup adaptation and
  ignores it.

- progress:

  Show sampler progress. \`FALSE\` fully silences it for every backend:
  the Turing \`Sampling ... ETA\` bar (via \`Turing.setprogress!\`), and
  Stan's per-iteration refresh, chain messages, and auto-opened progress
  window. (The MLE optimiser's progress is controlled separately, by
  \[optim_control()\]'s \`progress\`.)

- sigma_prior:

  Prior for the Gaussian/lognormal noise SD, as a Julia expression
  (default a half-Normal). Used by the Julia backend.

- phi_prior:

  Prior for the negative-binomial dispersion, as a Julia expression.
  Used by the Julia backend.

- sigma_prior_stan, phi_prior_stan:

  The same two priors written as \*\*Stan\*\* expressions (defaults
  \`"normal(0, 1)"\` and \`"gamma(2, 0.2)"\`). Used by the Stan backend,
  which cannot read the Julia-syntax versions above.

- adapt_delta:

  Stan NUTS target acceptance probability (default 0.8). Raise toward 1
  (e.g. 0.95) to reduce divergent transitions. Stan backend only.

## Value

A named list of Bayesian settings, passed as \`bayes =\` to
\[fitCompartmentalModel()\].

## See also

\[fitCompartmentalModel()\], \[solver_control()\], \[optim_control()\].

## Examples

``` r
bayes_control(chains = 4, iter = 2000, warmup = 1000)
#> $sampler
#> [1] "NUTS(0.65)"
#> 
#> $chains
#> [1] 4
#> 
#> $iter
#> [1] 2000
#> 
#> $warmup
#> [1] 1000
#> 
#> $seed
#> NULL
#> 
#> $init_from_optim
#> [1] TRUE
#> 
#> $progress
#> [1] TRUE
#> 
#> $sigma_prior
#> [1] "truncated(Normal(0, 1), 0, Inf)"
#> 
#> $phi_prior
#> [1] "Gamma(2, 5)"
#> 
#> $sigma_prior_stan
#> [1] "normal(0, 1)"
#> 
#> $phi_prior_stan
#> [1] "gamma(2, 0.2)"
#> 
#> $adapt_delta
#> [1] 0.8
#> 
bayes_control(progress = FALSE)              # silent sampling
#> $sampler
#> [1] "NUTS(0.65)"
#> 
#> $chains
#> [1] 4
#> 
#> $iter
#> [1] 2000
#> 
#> $warmup
#> [1] 1000
#> 
#> $seed
#> NULL
#> 
#> $init_from_optim
#> [1] TRUE
#> 
#> $progress
#> [1] FALSE
#> 
#> $sigma_prior
#> [1] "truncated(Normal(0, 1), 0, Inf)"
#> 
#> $phi_prior
#> [1] "Gamma(2, 5)"
#> 
#> $sigma_prior_stan
#> [1] "normal(0, 1)"
#> 
#> $phi_prior_stan
#> [1] "gamma(2, 0.2)"
#> 
#> $adapt_delta
#> [1] 0.8
#> 
bayes_control(chains = 2, iter = 1000, seed = 1)  # reproducible, lighter run
#> $sampler
#> [1] "NUTS(0.65)"
#> 
#> $chains
#> [1] 2
#> 
#> $iter
#> [1] 1000
#> 
#> $warmup
#> [1] 1000
#> 
#> $seed
#> [1] 1
#> 
#> $init_from_optim
#> [1] TRUE
#> 
#> $progress
#> [1] TRUE
#> 
#> $sigma_prior
#> [1] "truncated(Normal(0, 1), 0, Inf)"
#> 
#> $phi_prior
#> [1] "Gamma(2, 5)"
#> 
#> $sigma_prior_stan
#> [1] "normal(0, 1)"
#> 
#> $phi_prior_stan
#> [1] "gamma(2, 0.2)"
#> 
#> $adapt_delta
#> [1] 0.8
#> 
```
