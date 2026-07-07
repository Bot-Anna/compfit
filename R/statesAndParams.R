# Parse one entry column (States or Parameters) into a typed intermediate form.
# Both columns share the same grammar, so this is the single parser for it
# (statesAndParams used to carry two near-identical copies -- the source of the
# name/length prior bug, which the States half still latently had):
#   *name=value          fixed (a number, or an expression of other names)
#   name=[lo,hi]         fitted, Uniform prior / box search
#   name=Dist(args)[..]  fitted, a named prior (box via .mle_box for the search)
#   name=[lo,hi]|init    fitted with a known initial value
#   name=<expr>          NOT a prior -> a function of other names
# Returns a list of typed pieces; statesAndParams() maps them to its outputs.
.parse_entry_group <- function(vec) {
  vec <- vec[!is.na(vec) & vec != ""]
  vec <- gsub(" ", "", vec)

  fixed  <- grep("^\\*", vec, value = TRUE)
  fitted <- grep("^\\*", vec, value = TRUE, invert = TRUE)

  # Fixed: *name=value (value may be numeric or an expression of other names).
  fixed_final  <- setNames(sub("^\\*.*?=(.*)", "\\1", fixed),
                           sub("^\\*(.*?)=.*",  "\\1", fixed))
  fixed_nonnum <- is.na(suppressWarnings(as.numeric(fixed_final)))

  # Fitted split: known-initial (name=[lo,hi]|init) vs without.
  with_v    <- fitted[grepl("\\|", fitted)]
  without_v <- fitted[!grepl("\\|", fitted)]

  # Without-init PRIOR entries (box or named distribution) -> a finite L-BFGS-B
  # box via .mle_box. Names come from the SAME filtered set as the values, so a
  # distributional prior can never inflate the name vector past the values
  # (the bug that used to hit both States and Parameters).
  prior_v <- without_v[vapply(without_v, .is_prior_entry, logical(1))]
  w_names <- sub("=.*", "", prior_v)
  boxes   <- lapply(sub("^[^=]*=", "", prior_v), function(rhs) .mle_box(parsePrior(rhs)))
  without_lower <- setNames(vapply(boxes, function(b) b[1], numeric(1)), w_names)
  without_upper <- setNames(vapply(boxes, function(b) b[2], numeric(1)), w_names)

  # With-init entries: name=[lo,hi]|init (box only). The suffix after '|' is
  # either a number (an explicit start) or the keyword 'random' / 'rand' (draw
  # the start uniformly from the box at fit time, seeded for reproducibility).
  wi_names   <- sub("=.*", "", with_v)
  # NA bounds are the intended sentinel for a non-box '|' entry (caught below /
  # upstream by validate), so the coercion is deliberately quiet.
  with_lower <- setNames(suppressWarnings(as.numeric(sub(".*=\\[([^]]+),.*", "\\1", with_v))), wi_names)
  with_upper <- setNames(suppressWarnings(as.numeric(sub(".*\\[.*,(\\d+(\\.\\d+)?)\\].*", "\\1", with_v))), wi_names)
  with_suffix <- sub(".*\\|", "", with_v)
  with_random <- setNames(grepl("^(random|rand)$", with_suffix, ignore.case = TRUE), wi_names)
  # A random start is defined by the box it draws from, so it needs a box prior.
  bad_rnd <- with_random & (is.na(with_lower) | is.na(with_upper))
  if (any(bad_rnd))
    stop("A '|random' start needs a box prior, '[lo,hi]|random': ",
         paste(with_v[bad_rnd], collapse = ", "), call. = FALSE)
  with_init  <- setNames(suppressWarnings(as.numeric(with_suffix)), wi_names)
  # Random entries get the box midpoint as a placeholder init; the actual
  # uniform draw happens at fit time (in .fit_optim, after set.seed) so it is
  # reproducible. The placeholder keeps a finite value for the normalisation.
  with_init[with_random] <- (with_lower[with_random] + with_upper[with_random]) / 2

  # Fitted entries that are NOT priors -> expressions of other names (functions).
  # These are non-'*' entries, so name = before '=', value = after (the old code
  # used a '^\\*' regex that never matched here, mangling both).
  fn_v      <- fitted[!vapply(fitted, .is_prior_entry, logical(1))]
  functions <- setNames(sub("^[^=]*=", "", fn_v), sub("=.*", "", fn_v))

  list(fixed = fixed_final, fixed_nonnumeric = fixed_nonnum,
       without_names = w_names, without_lower = without_lower, without_upper = without_upper,
       with_names = wi_names, with_lower = with_lower, with_upper = with_upper, with_init = with_init,
       with_random = with_random, functions = functions)
}

statesAndParams <- function(modelParams) {
  comp_names <- .compartments(modelParams)   # canonical compartment order (names)
  S <- .parse_entry_group(modelParams$States)
  P <- .parse_entry_group(modelParams$Parameters)

  ## ---- States ----
  states_fixed_final <- S$fixed
  non_numeric_fixed  <- S$fixed_nonnumeric
  # Expression-valued fixed states (parameter-dependent initial states) become
  # closures of `parms`: substitute each free variable v with parms[[v]].
  exprs <- states_fixed_final[non_numeric_fixed]
  list_states_functions <- lapply(exprs, function(txt) {
    e <- parse(text = txt)[[1]]
    replaced <- e
    for (v in all.vars(e))
      replaced <- do.call("substitute",
                          list(replaced, setNames(list(substitute(parms[[v]])), v)))
    eval(call("function", as.pairlist(list(parms = NULL)), replaced))
  })
  names(list_states_functions) <- names(exprs)

  lower_guesses_states_without <- S$without_lower
  upper_guesses_states_without <- S$without_upper
  lower_guesses_states_with    <- S$with_lower
  upper_guesses_states_with    <- S$with_upper
  initial_states_with          <- S$with_init
  states_fitted_limits_2 <- c(setNames(S$with_names,    S$with_names),
                              setNames(S$without_names, S$without_names))
  lower_guesses_states <- c(lower_guesses_states_with, lower_guesses_states_without)
  upper_guesses_states <- c(upper_guesses_states_with, upper_guesses_states_without)
  states_functions     <- gsub("\\btime\\b", "t", S$functions)   # Julia: time -> t

  ## ---- Parameters ----
  params_fixed_final <- P$fixed
  lower_guesses_params_without <- P$without_lower
  upper_guesses_params_without <- P$without_upper
  lower_guesses_params_with    <- P$with_lower
  upper_guesses_params_with    <- P$with_upper
  initial_params_with          <- P$with_init
  params_fitted_limits_2 <- c(setNames(P$with_names,    P$with_names),
                              setNames(P$without_names, P$without_names))
  lower_guesses_params <- c(lower_guesses_params_with, lower_guesses_params_without)
  upper_guesses_params <- c(upper_guesses_params_with, upper_guesses_params_without)
  params_functions     <- gsub("\\btime\\b", "t", P$functions)

  ## ---- Initial guesses / hypercube / not-Sobol normalisation ----
  initial_guesses <- c(initial_states_with,
                       midpoint(lower_guesses_states_without, upper_guesses_states_without),
                       initial_params_with,
                       midpoint(lower_guesses_params_without, upper_guesses_params_without))
  dimension_of_hypercube <- length(lower_guesses_states_without) +
                            length(lower_guesses_params_without)
  initial_guesses_not_sobol_states_normalised <- mapply(
    normalise, initial_states_with, lower_guesses_states_with, upper_guesses_states_with)
  initial_guesses_not_sobol_params_normalised <- mapply(
    normalise, initial_params_with, lower_guesses_params_with, upper_guesses_params_with)

  ## ---- Random-start mask ----
  # Which fitted quantities requested a '|random' start, aligned to the fitted
  # order (with-init entries first, then without-init, matching lower_*/upper_*).
  # Without-init entries never carry '|', so they are always FALSE.
  random_states <- c(S$with_random,
                     setNames(rep(FALSE, length(S$without_names)), S$without_names))
  random_params <- c(P$with_random,
                     setNames(rep(FALSE, length(P$without_names)), P$without_names))

  list(
    states_fitted        = states_fitted_limits_2,
    upper_states         = upper_guesses_states,
    lower_states         = lower_guesses_states,
    states_fixed         = states_fixed_final,
    states_functions     = states_functions,
    params_fitted        = params_fitted_limits_2,
    upper_params         = upper_guesses_params,
    lower_params         = lower_guesses_params,
    params_fixed         = params_fixed_final,
    params_functions     = params_functions,
    initial_guesses      = initial_guesses,
    hypercube_dim        = dimension_of_hypercube,
    not_sobol_states     = initial_guesses_not_sobol_states_normalised,
    not_sobol_params     = initial_guesses_not_sobol_params_normalised,
    list_states_functions = list_states_functions,
    non_numeric_fixed    = non_numeric_fixed,
    lower_states_without = lower_guesses_states_without,
    upper_states_without = upper_guesses_states_without,
    lower_params_without = lower_guesses_params_without,
    upper_params_without = upper_guesses_params_without,
    random_states        = random_states,
    random_params        = random_params,
    comp_names           = comp_names
  )
}
