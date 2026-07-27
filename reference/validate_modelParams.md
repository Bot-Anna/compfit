# Validate a modelParams sheet

Check a model-parameter sheet for common entry errors before it is
built. On any problem it raises a single error listing every issue
found, each naming the column/cell and what is wrong. Checks:

- fixed `States`/`Parameters` values and `Linear`/
  `Quadratic`/`Constant`/`Functions`/`Conditions` expressions reference
  only declared symbols (parameters, `_0` aliases, the declared states –
  `X1..Xn` or named, e.g. `S`/`I`/`R` –, functions, `time`), so a typo
  like `*2*-bta` is caught;

- box priors `[lo,hi]` are two numbers with `lo < hi`;

- distribution priors have numeric arguments;

- no parameter/state/function name collides with a Julia keyword or a
  codegen variable/function (`t`/`p`/`X`/`du`/`dX`/ `parms`, `N1..`,
  `total_pop`, `f<ij>`, `g<ij>`, `cst<i>`, `secOrd_i_j`, or another
  quantity's `_0` alias);

- `Quadratic` cells are `*goto*coeff` with a target in `1..n`;

- `Others` has numeric `startpoint`/`endpoint`/`partition`.

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
