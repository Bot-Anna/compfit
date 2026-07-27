# Export a fit as self-contained source code

Emits portable R source that reproduces the model, loss, and/or fit from
a \`"compartmentalFit"\` (or model) object, inlining the captured data
so the script runs in a fresh session (given a Julia bridge).

## Usage

``` r
extract_code(
  x,
  file = NULL,
  what = c("all", "model", "loss", "fit"),
  inline = c("readable", "fidelity")
)
```

## Arguments

- x:

  A \`"compartmentalFit"\` or \`"compartmentalModel"\` object.

- file:

  Optional output path; if \`NULL\`, returns the code as a string.

- what:

  Which part(s) to emit: all, model, loss, or fit.

- inline:

  How captured data objects (bounds, data matrices, masks,
  \`partition\`, \`init_norm\`, ...) are written into the script.
  \`"readable"\` (default) emits them as human-readable, re-parseable
  \`deparse()\` code (names, dims, NAs preserved; doubles to 17
  significant digits, which round-trips exactly for typical values but
  can differ by ~1 ULP on adversarial doubles, as R's decimal parser is
  not always correctly rounded). \`"fidelity"\` emits them as opaque
  \`unserialize(as.raw(...))\` blobs — a guaranteed byte-exact
  round-trip, for when the last bit matters.

## Value

The generated code as a character string, invisibly (also written to
\`file\` when supplied).

## Examples

``` r
if (FALSE) { # \dontrun{
extract_code(fit, file = file.path(tempdir(), "reproduce.R"), what = "all")
# opaque-but-byte-exact data blocks:
extract_code(fit, file = file.path(tempdir(), "exact.R"), inline = "fidelity")
} # }
```
