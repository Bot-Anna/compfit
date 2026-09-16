# Replay a model's derived quantities on a solved trajectory

Internal, but EXPORTED deliberately: \[extract_code()\] emits a
self-contained script whose loss body calls \`.derived_columns()\`
unqualified after a plain \`library(compfit)\`, so it must be reachable
from outside the namespace. The call there sits inside a
\`tryCatch(error = NULL)\` that degrades to states-and-parameters only,
so un-exporting this would not raise an error in the extracted script –
it would silently drop every derived column from the loss. Not part of
the stable user-facing API.

## Usage

``` r
.derived_columns(sir_out, parms, time, spec)
```

## Arguments

- sir_out:

  Solved trajectory: one column per state, over the time grid.

- parms:

  Named numeric vector of parameters.

- time:

  Numeric time grid.

- spec:

  The \`derived_spec\` built by \`compartmentalFunction()\`.

## Value

A data frame of derived columns, or \`NULL\` if there is nothing to add.
