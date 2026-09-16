# Central point estimate of a fit, whichever method produced it

Returns \`list(initial_state, parms)\` on the natural scale for ANY
successful fit: the optimiser's point estimate for an MLE-type fit, or
the posterior means for \`method = "bayes"\` (where \`fit\$point\` is
\`NULL\`, since sampling produces a distribution rather than a point).
This is the vector \`plot_fit()\` draws its central trajectory at, so it
is also what to hand to \[solve_and_evaluate()\] to reproduce or extend
that trajectory.

## Usage

``` r
get_central_point(fit)
```

## Arguments

- fit:

  A \`"compartmentalFit"\` object.

## Value

A list with \`initial_state\` and \`parms\`, both named numeric vectors.

## Details

For a Bayesian fit the posterior mean is a summary, not a sampled draw:
it can sit in a low-probability region when the posterior is skewed or
multimodal, and for a strongly non-linear model the trajectory at the
mean parameters is not the mean trajectory. Use \[posterior_draws()\] to
work draw by draw.

## See also

\[get_point()\] (MLE-type fits only), \[posterior_means()\],
\[solve_and_evaluate()\].

## Examples

``` r
if (FALSE) { # \dontrun{
p <- get_central_point(fit)
ev <- solve_and_evaluate(fit, p$initial_state, p$parms)
} # }
```
