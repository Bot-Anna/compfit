# Hypercube pre-search settings

Settings for the quasi-Monte Carlo (Sobol) hypercube pre-search
(\`method = "hypercube"\`), which evaluates the loss on a
low-discrepancy Sobol design over the normalised \\0,1\\ box and keeps
the best point. Deterministic and derivative-free – useful as a
stand-alone coarse fit or to seed a warm start for the local optimiser.

## Usage

``` r
hypercube_control(n = 10000)
```

## Arguments

- n:

  Number of Sobol sample points to evaluate (default 10000). More points
  give finer coverage of the parameter box at linear cost in solves. n
  should approximately be of magnitude 10^m where m is the number of
  unknown parameters and initial states.

## Value

A named list with the sample count, passed as \`hypercube =\` to
\[fitCompartmentalModel()\].

## See also

\[fitCompartmentalModel()\], \[optim_control()\].

## Examples

``` r
hypercube_control(n = 5000)
#> $n
#> [1] 5000
#> 
```
