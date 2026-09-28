# Move a fitted model's start year earlier, keeping its parameter estimates

Rewrites a fit's model sheet so the time grid begins at an earlier
\`startpoint\`, with every fitted parameter frozen at its estimate
(\`\*name=value\`, via \[fill_params()\]) and the data padded with blank
columns for the new years. The initial state then refers to the NEW
start year, so it is almost always wrong as inherited and should be
respecified through \`states\`.

## Usage

``` r
reanchor_scenario(
  fit,
  startpoint,
  states = NULL,
  modelParams = fit$model$modelParams,
  dataCombined = fit$data$data_combined,
  summary = c("median", "mean"),
  digits = NULL
)
```

## Arguments

- fit:

  A successful \`"compartmentalFit"\` object.

- startpoint:

  New (earlier) start year. Must be a whole number strictly before the
  fit's current \`startpoint\`; to go the other way see the \`endpoint\`
  argument of \[plot_fit()\] and \[solve_and_evaluate()\].

- states:

  Optional named character vector respecifying the initial state at the
  new start year, e.g. \`c(X1 = "\[900,1000\]", X2 = "\[0,50\]", X3 =
  "0")\`. Entries follow the sheet's own grammar: a box or distribution
  is FITTED, a bare number or an expression is fixed. Names must be
  declared compartments. \`NULL\` (the default) inherits the filled
  state and warns, since that describes the original start year.

- modelParams:

  Model sheet to rewrite (defaults to the fit's own).

- dataCombined:

  Observed data to pad (defaults to the fit's own).

- summary:

  Which posterior summary to freeze for a Bayesian fit, \`"median"\`
  (default) or \`"mean"\`; passed to \[fill_params()\].

- digits:

  Optional rounding for the frozen values.

## Value

A list with \`modelParams\` and \`dataCombined\`, ready for
\[fitCompartmentalModel()\] or \[simulate_model()\].

## Details

This is re-anchoring, not backward integration: the model is solved
forward from the earlier year, so the earlier initial state is a
modelling assumption you supply (or estimate), not something recovered
from the fit. Nothing is solved here – the returned sheet must be
rebuilt, because the generated ODE carries \`startpoint\` as a literal.

## See also

\[fill_params()\], \[simulate_model()\], \[fitCompartmentalModel()\],
\[plot_fit()\] (its \`endpoint\` argument projects forward instead).

## Examples

``` r
if (FALSE) { # \dontrun{
re <- reanchor_scenario(fit, startpoint = 2010,
        states = c(X1 = "[900,1000]", X2 = "[0,50]", X3 = "0"))
early <- fitCompartmentalModel(re$modelParams, re$dataCombined,
           solver = solver_control(backend = "r"))
} # }
```
