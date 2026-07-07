# numberOfComps returns a list with
# - the compartment-index columns (any column whose name starts with "_",
#   e.g. _Level1, _Level2, ...)
# - the compartment indices they hold (by number), mapped to canonical position
# - the canonical compartment names (States-column order)
# - total number of compartments
numberOfComps <- function(modelParams) {
  compartment_cols <- grep("^_", names(modelParams),
                       ignore.case = TRUE, value = TRUE)
  # Canonical registry: the States column defines the compartments AND their
  # order (compartment i = i-th States entry). _Level cells may name a
  # compartment either by number (1..n) or by that name, so both resolve here.
  comp_names <- .compartments(modelParams)
  compartment_values <- unlist(
    lapply(modelParams[compartment_cols], function(col) {
      col[col != "" & !is.na(col)]
    }),
    use.names = FALSE
  )
  # Map each _Level cell (a number OR a compartment name) to its canonical
  # index via the registry, so downstream membership tests (i %in% vec) and
  # X[j] subscripts are always plain integers -- "10" beats "9" numerically,
  # and "S" resolves to whatever position S holds in the States column.
  compartment_values <- .comp_index(compartment_values, comp_names)
  compartment_values <- compartment_values[!is.na(compartment_values)]
  # Prefer the registry for the count (States is the source of truth); fall
  # back to the largest _Level index for a legacy sheet without a States column.
  number_of_comps <- if (length(comp_names)) length(comp_names)
                     else max(compartment_values)

  return(list(compartment_cols = compartment_cols,
              compartment_values = compartment_values,
              comp_names = comp_names,
              number_of_comps = number_of_comps))
}