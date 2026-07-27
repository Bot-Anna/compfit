# Initialise the Julia bridge

Starts the Julia session (via diffeqr + JuliaCall) and, on first use,
provisions the required Julia packages (SciMLBase, OrdinaryDiffEq,
Turing, Distributions) from a *pinned* environment shipped with the
package (`inst/julia/Project.toml`), copied to a writable per-user cache
and instantiated so every user gets the same tested versions. If that
fails it falls back to installing the packages into the default Julia
environment. Loads OrdinaryDiffEq. Safe to call repeatedly: the session
bind and the one-time provisioning are each guarded.

## Usage

``` r
setup_julia()
```

## Value

TRUE, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
setup_julia()  # starts Julia; installs the Julia packages on first call
} # }
```
