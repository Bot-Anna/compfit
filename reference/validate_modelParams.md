# Validate a modelParams sheet

Check a model-parameter sheet for common entry errors before it is
built. On any problem it raises a single error listing every issue
found, each naming the column/cell and what is wrong. Checks:

- fixed `States`/`Parameters` values and `Linear`/
  `Quadratic`/`Constant`/`Functions`/`Conditions`/
  `Mixing_<level>`/`Pool_<level>` expressions reference only declared
  symbols (parameters, `_0` aliases, the declared states – `X1..Xn` or
  named, e.g. `S`/`I`/`R`, which are rewritten to their state slot –,
  functions, `time`, `N<level>`), so a typo like `*2*-bta` is caught; a
  parameter/function may not share a State name;

- box priors `[lo,hi]` are two numbers with `lo < hi`;

- distribution priors have numeric arguments;

- no parameter/state/function name collides with a Julia keyword or a
  codegen variable (`t`/`p`/`X`/`du`/`dX`/ `parms`, `N1..`, `Nw1..`,
  `total_pop`, `secOrd_i_j`, or another quantity's `_0` alias);

- every declared name (Parameters/Functions/States) is a legal
  identifier – letters, digits, `_`, `.` – so operator characters (`^`,
  `()`, `,`, `-`, `/`, ...) that would be parsed as operators rather
  than a name are rejected;

- `Quadratic` cells are `*goto*coeff` with a target given as a
  compartment index (`1..n`) or a State name (e.g. `*I*`);

- each `Mixing_<level>` / `Pool_<level>` column has a matching
  `_<level>` index column (`Mixing` holds non-negative per-compartment
  weights; `Pool` holds a single pool expression, and the two are
  mutually exclusive for a level);

- `Others` has numeric `startpoint`/`endpoint`/`partition`;

- the `Functions` column can be put in a working evaluation order. Row
  order itself does not matter – a Function may be written above the
  entries it references, and the builder sorts the column into
  dependency order – but a name **defined more than once**, a
  **self-reference** (`x <- x + 1`; a Function is a definition, not an
  update), and a **reference cycle** have no valid order and are
  rejected, the cycle reported as the loop it forms.

## Usage

``` r
validate_modelParams(modelParams)
```

## Arguments

- modelParams:

  The model-parameter data frame (as read by
  [`load_scenario()`](load_scenario.md)).

## Value

Invisibly `TRUE` when valid; otherwise stops with the collected
problems.

## Examples

``` r
mini <- system.file("extdata", "minimal", package = "compfit")
sc <- load_scenario(mini, combined_file = "dataCombined.csv",
                    dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
validate_modelParams(sc$modelParams)   # TRUE (the fixture is valid)
```
