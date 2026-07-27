# Optimiser settings

Bundle optimiser settings for the maximum-likelihood methods. Two
optimisers share this control: the deterministic quasi-Newton
\*\*L-BFGS-B\*\* (\`method = "lbfgsb"\`, gradient by finite differences)
and the stochastic \*\*DEoptim\*\* differential-evolution search
(\`method = "deoptim"\`, a global method that is more robust to
multi-modal / rugged objectives but slower). The \`maxit\`/\`factr\`
options apply to L-BFGS-B; \`itermax\`/\`NP_mult\`/\`CR\`/\`F\`/
\`trace\`/\`seed\` apply to DEoptim; \`progress\` applies to both.
Settings for the other optimiser are simply ignored, so one control
object works for any MLE method.

## Usage

``` r
optim_control(
  opt_method = "L-BFGS-B",
  maxit = 1000,
  factr = 1e+07,
  itermax = 1000,
  NP_mult = 10,
  CR = 0.9,
  F = 0.8,
  trace = 50,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- opt_method:

  The \`optim()\` method used on the \`lbfgsb\` path. Normally
  \`"L-BFGS-B"\` (box-constrained); \`"Brent"\` is the 1-D alternative.
  The search runs over a normalised \\0,1\\ box, so the method must
  support box bounds.

- maxit:

  L-BFGS-B: maximum iterations before it stops (default 1000).

- factr:

  L-BFGS-B: convergence tolerance as a multiple of machine epsilon;
  smaller = stricter/slower (default 1e7, i.e. ~1e-8 relative).

- itermax:

  DEoptim: maximum number of generations (default 1000).

- NP_mult:

  DEoptim: population size as a multiple of the number of fitted
  quantities (\`NP = NP_mult \* d\`). The DEoptim authors recommend \>=
  10 (the default); lower is faster but explores less.

- CR:

  DEoptim: crossover probability in \\0,1\\ (default 0.9).

- F:

  DEoptim: differential weighting (mutation) factor, typically in
  \\0,2\\ (default 0.8).

- trace:

  DEoptim: print the best value every \`trace\` generations (default
  50). Ignored when \`progress = FALSE\`.

- seed:

  Optional RNG seed for the stochastic \`deoptim\` method, for
  reproducible fits. \`lbfgsb\`/\`hypercube\` are already deterministic
  and ignore it – except that a \`\|random\` parameter start (see the
  modelParams grammar in \`vignette("compfit")\`) draws its L-BFGS-B
  start from this seed, so setting it makes that random warm start
  reproducible too.

- progress:

  Show optimisation progress. \`TRUE\` (default) prints the running best
  objective each time the loss improves and (for DEoptim) its
  per-generation \`trace\`. \`FALSE\` silences BOTH for a quiet fit –
  the best solution is still tracked and recovered. (The Bayesian
  sampler's progress is controlled separately, by \[bayes_control()\]'s
  \`progress\`.)

## Value

A named list of optimiser settings, passed as \`control =\` to
\[fitCompartmentalModel()\].

## See also

\[fitCompartmentalModel()\], \[bayes_control()\],
\[hypercube_control()\].

## Examples

``` r
optim_control(maxit = 500)                       # L-BFGS-B, 500 iterations
#> $opt_method
#> [1] "L-BFGS-B"
#> 
#> $seed
#> NULL
#> 
#> $maxit
#> [1] 500
#> 
#> $factr
#> [1] 1e+07
#> 
#> $itermax
#> [1] 1000
#> 
#> $NP_mult
#> [1] 10
#> 
#> $CR
#> [1] 0.9
#> 
#> $F
#> [1] 0.8
#> 
#> $trace
#> [1] 50
#> 
#> $progress
#> [1] TRUE
#> 
optim_control(opt_method = "DEoptim", itermax = 200)  # DE, 200 generations
#> $opt_method
#> [1] "DEoptim"
#> 
#> $seed
#> NULL
#> 
#> $maxit
#> [1] 1000
#> 
#> $factr
#> [1] 1e+07
#> 
#> $itermax
#> [1] 200
#> 
#> $NP_mult
#> [1] 10
#> 
#> $CR
#> [1] 0.9
#> 
#> $F
#> [1] 0.8
#> 
#> $trace
#> [1] 50
#> 
#> $progress
#> [1] TRUE
#> 
optim_control(progress = FALSE)                  # quiet fit (no console output)
#> $opt_method
#> [1] "L-BFGS-B"
#> 
#> $seed
#> NULL
#> 
#> $maxit
#> [1] 1000
#> 
#> $factr
#> [1] 1e+07
#> 
#> $itermax
#> [1] 1000
#> 
#> $NP_mult
#> [1] 10
#> 
#> $CR
#> [1] 0.9
#> 
#> $F
#> [1] 0.8
#> 
#> $trace
#> [1] 50
#> 
#> $progress
#> [1] FALSE
#> 
optim_control(seed = 1)                          # reproducible DEoptim
#> $opt_method
#> [1] "L-BFGS-B"
#> 
#> $seed
#> [1] 1
#> 
#> $maxit
#> [1] 1000
#> 
#> $factr
#> [1] 1e+07
#> 
#> $itermax
#> [1] 1000
#> 
#> $NP_mult
#> [1] 10
#> 
#> $CR
#> [1] 0.9
#> 
#> $F
#> [1] 0.8
#> 
#> $trace
#> [1] 50
#> 
#> $progress
#> [1] TRUE
#> 
```
