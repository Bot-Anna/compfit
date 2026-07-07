# ============================================================
# priorSpec.R
# Parallel prior-specification builder (Bayesian path only).
#
# This reads modelParams INDEPENDENTLY of statesAndParams() and produces a
# structured prior spec per estimated quantity, in the SAME ORDER that the
# Julia ODE expects its parameter vector `p` and initial-state vector `X`.
# It does not modify or depend on statesAndParams()'s return value; it is a
# parallel structure, so the MLE path is entirely untouched.
#
# Ordering contract (must match generateExpressions.R / the Julia ODE):
#   p = [ fitted params..., fixed params... ]   (params_fitted first)
#   X = [ X1_0, X2_0, ... ]                      (compartment order)
# For PRIORS we only need the *estimated* quantities (fitted params + fitted
# states). Fixed values are returned separately so the Julia model can splice
# them into p / X at the right positions.
#
#   buildPriorSpec(modelParams)
#     -> list(
#          params = list(<name> = <parsePrior result>, ...),  # estimated, fitted-param order
#          states = list(<name> = <parsePrior result>, ...),  # estimated, fitted-state order
#          fixed_params = named numeric/char (as in the sheet),
#          fixed_states = named numeric/char,
#          order = list(params_fitted = chr, states_fitted = chr,
#                       params_fixed = chr, states_fixed = chr)
#        )
# ============================================================

buildPriorSpec <- function(modelParams) {
  
  ## ---- Helper: split a sheet column into fixed vs fitted entries ----
  clean_col <- function(col) {
    v <- col[col != "" & !is.na(col)]
    gsub(" ", "", v)
  }
  
  ## ---- Parameters ----
  params_vector <- clean_col(modelParams$Parameters)
  params_fixed_entries  <- grep("^\\*", params_vector, value = TRUE)
  params_fitted_entries <- grep("^\\*", params_vector, value = TRUE, invert = TRUE)
  
  # Fixed parameters: name=value (value may be a number or an expression string,
  # exactly as statesAndParams handles them). We keep them as-is.
  pf_names  <- sub("^\\*(.*?)=.*", "\\1", params_fixed_entries)
  pf_values <- sub("^\\*.*?=(.*)", "\\1", params_fixed_entries)
  fixed_params <- setNames(pf_values, pf_names)
  
  # Fitted parameters: name=<prior-rhs>. parsePrior interprets the rhs.
  # Only those that actually have a "[" or a distribution are estimated; a fitted
  # entry that is itself a function of others (no "=[", no Dist()) is NOT a prior
  # target and is left for the existing function-handling machinery. We detect a
  # prior target as "has =[ or =Dist(".
  is_prior_target <- function(entry) {
    rhs <- sub("^[^=]*=", "", entry)
    grepl("^\\[", rhs) || grepl("^[A-Za-z]+\\(", rhs)
  }
  pt <- params_fitted_entries[vapply(params_fitted_entries, is_prior_target, logical(1))]
  pt_names <- sub("=.*", "", pt)
  params_spec <- setNames(
    lapply(pt, function(entry) parsePrior(sub("^[^=]*=", "", entry))),
    pt_names
  )
  
  ## ---- States ----
  states_vector <- clean_col(modelParams$States)
  states_fixed_entries  <- grep("^\\*", states_vector, value = TRUE)
  states_fitted_entries <- grep("^\\*", states_vector, value = TRUE, invert = TRUE)
  
  sf_names  <- sub("^\\*(.*?)=.*", "\\1", states_fixed_entries)
  sf_values <- sub("^\\*.*?=(.*)", "\\1", states_fixed_entries)
  fixed_states <- setNames(sf_values, sf_names)
  
  st <- states_fitted_entries[vapply(states_fitted_entries, is_prior_target, logical(1))]
  st_names <- sub("=.*", "", st)
  states_spec <- setNames(
    lapply(st, function(entry) parsePrior(sub("^[^=]*=", "", entry))),
    st_names
  )
  
  list(
    params       = params_spec,
    states       = states_spec,
    fixed_params = fixed_params,
    fixed_states = fixed_states,
    order = list(
      params_fitted = pt_names,
      states_fitted = st_names,
      params_fixed  = pf_names,
      states_fixed  = sf_names
    )
  )
}