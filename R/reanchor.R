# ============================================================
# reanchor.R
# reanchor_scenario() -- move a fitted model's startpoint EARLIER and keep its
# parameter estimates, so the model can be run over years before the data.
#
# This is the well-posed way to "project backwards". Integrating the ODE in
# reverse from the fitted initial state is not: forward epidemic dynamics are
# contracting, so running them backwards is expanding, and the solution leaves
# the feasible region (negative compartments, susceptibles above the population).
# Instead we re-anchor: keep the fitted parameters, move the start year back, and
# treat the initial state at the NEW start year as fresh unknowns -- fixed by hand
# or fitted against whatever early data exists. Integration always runs forward.
#
# The result is a sheet + data pair, not a fit, because the model must be
# regenerated: `startpoint` is emitted into the generated ODE as a literal (unlike
# `endpoint`, which is why plot_fit(endpoint=) can reuse a registered model). Hand
# the pair to fitCompartmentalModel() to estimate the earlier initial state, or to
# simulate_model() if you fixed every quantity.
# ============================================================

#' Move a fitted model's start year earlier, keeping its parameter estimates
#'
#' Rewrites a fit's model sheet so the time grid begins at an earlier
#' `startpoint`, with every fitted parameter frozen at its estimate
#' (`*name=value`, via [fill_params()]) and the data padded with blank columns for
#' the new years. The initial state then refers to the NEW start year, so it is
#' almost always wrong as inherited and should be respecified through `states`.
#'
#' This is re-anchoring, not backward integration: the model is solved forward
#' from the earlier year, so the earlier initial state is a modelling assumption
#' you supply (or estimate), not something recovered from the fit. Nothing is
#' solved here -- the returned sheet must be rebuilt, because the generated ODE
#' carries `startpoint` as a literal.
#'
#' @param fit A successful `"compartmentalFit"` object.
#' @param startpoint New (earlier) start year. Must be a whole number strictly
#'   before the fit's current `startpoint`; to go the other way see the
#'   `endpoint` argument of [plot_fit()] and [solve_and_evaluate()].
#' @param states Optional named character vector respecifying the initial state
#'   at the new start year, e.g.
#'   `c(X1 = "[900,1000]", X2 = "[0,50]", X3 = "0")`. Entries follow the sheet's
#'   own grammar: a box or distribution is FITTED, a bare number or an expression
#'   is fixed. Names must be declared compartments. `NULL` (the default) inherits
#'   the filled state and warns, since that describes the original start year.
#' @param modelParams Model sheet to rewrite (defaults to the fit's own).
#' @param dataCombined Observed data to pad (defaults to the fit's own).
#' @param summary Which posterior summary to freeze for a Bayesian fit,
#'   `"median"` (default) or `"mean"`; passed to [fill_params()].
#' @param digits Optional rounding for the frozen values.
#' @return A list with `modelParams` and `dataCombined`, ready for
#'   [fitCompartmentalModel()] or [simulate_model()].
#' @seealso [fill_params()], [simulate_model()], [fitCompartmentalModel()],
#'   [plot_fit()] (its `endpoint` argument projects forward instead).
#' @examples
#' \dontrun{
#' re <- reanchor_scenario(fit, startpoint = 2010,
#'         states = c(X1 = "[900,1000]", X2 = "[0,50]", X3 = "0"))
#' early <- fitCompartmentalModel(re$modelParams, re$dataCombined,
#'            solver = solver_control(backend = "r"))
#' }
#' @export
reanchor_scenario <- function(fit, startpoint, states = NULL,
                              modelParams = fit$model$modelParams,
                              dataCombined = fit$data$data_combined,
                              summary = c("median", "mean"), digits = NULL) {
  .is_fit(fit)
  summary <- match.arg(summary)
  if (is.null(modelParams))
    stop("reanchor_scenario(): no model sheet on the fit; pass modelParams =.",
         call. = FALSE)

  tg <- .time_grid(modelParams)
  if (!is.numeric(startpoint) || length(startpoint) != 1L || !is.finite(startpoint))
    stop("'startpoint' must be a single finite year (got: ",
         paste(format(startpoint), collapse = ", "), ").", call. = FALSE)
  if (startpoint != round(startpoint))
    stop(sprintf("'startpoint' must be a whole year (got: %g). The time axis is annual.",
                 startpoint), call. = FALSE)
  startpoint <- as.numeric(round(startpoint))
  if (startpoint >= tg$startpoint)
    stop(sprintf(paste0(
      "'startpoint' (%g) must be EARLIER than the fitted startpoint (%g). ",
      "reanchor_scenario() extends a model backwards; to run past the fitted ",
      "endpoint use plot_fit(endpoint =) or solve_and_evaluate(endpoint =)."),
      startpoint, tg$startpoint), call. = FALSE)

  # Freeze the estimates, then move the start year.
  mp <- fill_params(fit, modelParams, digits = digits, summary = summary)
  mp <- .ra_set_others(mp, "startpoint", startpoint)

  # Respecify the initial state at the new start year.
  comps <- .compartments(mp)
  if (!is.null(states)) {
    if (is.null(names(states)) || any(!nzchar(names(states))))
      stop("'states' must be a NAMED vector, e.g. c(X1 = \"[900,1000]\").", call. = FALSE)
    bad <- setdiff(names(states), comps)
    if (length(bad))
      stop(sprintf("'states' names %s are not declared compartments (%s).",
                   paste(sQuote(bad), collapse = ", "),
                   paste(comps, collapse = ", ")), call. = FALSE)
    mp <- .ra_set_states(mp, states, comps)
  } else {
    warning(sprintf(paste0(
      "The initial state was inherited from the fit, so it still describes %g, ",
      "not %g. Respecify it with states = (e.g. c(%s = \"[lo,hi]\")) unless that ",
      "is deliberate."), tg$startpoint, startpoint, comps[1]), call. = FALSE)
  }

  list(modelParams = mp,
       dataCombined = .ra_pad_data(dataCombined, startpoint, tg$startpoint))
}


# Replace (or add) a `key=value` entry in the Others column, preserving the rest.
.ra_set_others <- function(mp, key, value) {
  others <- mp$Others
  hit <- which(grepl(sprintf("^\\s*%s\\s*=", key), others))
  entry <- sprintf("%s=%s", key, format(value, scientific = FALSE))
  if (length(hit)) {
    others[hit[1]] <- entry
  } else {
    slot <- which(is.na(others) | !nzchar(others))
    if (!length(slot))
      stop(sprintf("No free row in 'Others' to add '%s'; add a blank row.", entry),
           call. = FALSE)
    others[slot[1]] <- entry
  }
  mp$Others <- others
  mp
}


# Rewrite the States cells named in `states`. A box/distribution entry is FITTED
# and must NOT carry the '*' fixed marker; a bare number or an expression is
# fixed, and an expression needs the marker.
.ra_set_states <- function(mp, states, comps) {
  cells <- mp$States
  for (nm in names(states)) {
    val <- trimws(as.character(states[[nm]]))
    is_box  <- grepl("^\\[", val) || grepl("^[A-Za-z]+\\(", val)
    is_num  <- !is.na(suppressWarnings(as.numeric(val)))
    cell <- if (is_box || is_num) sprintf("%s=%s", nm, val)
            else sprintf("*%s=%s", nm, val)
    row <- which(vapply(cells, function(c0) {
      if (is.na(c0) || !nzchar(c0)) return(FALSE)
      identical(sub("=.*", "", sub("^\\*", "", gsub(" ", "", c0))), nm)
    }, logical(1)))
    if (!length(row))
      stop(sprintf("No States cell found for compartment '%s'.", nm), call. = FALSE)
    cells[row[1]] <- cell
  }
  mp$States <- cells
  mp
}


# Pad the data with BLANK columns for the new years. Blank cells are 'missing' in
# the data-cell grammar, so the earlier years contribute nothing to the loss --
# no data is invented, and .prepare_data's column-count-vs-horizon guard passes.
.ra_pad_data <- function(dataCombined, new_start, old_start) {
  if (is.null(dataCombined)) return(NULL)
  dc <- as.data.frame(dataCombined, stringsAsFactors = FALSE, check.names = FALSE)
  meta <- intersect(.DATA_META_COLS, names(dc))
  years <- as.character(seq(new_start, old_start - 1))
  years <- setdiff(years, names(dc))               # never clobber an existing column
  if (!length(years)) return(dc)
  pad <- setNames(replicate(length(years), rep(NA_character_, nrow(dc)),
                            simplify = FALSE), years)
  cbind(dc[meta], as.data.frame(pad, check.names = FALSE),
        dc[setdiff(names(dc), meta)], stringsAsFactors = FALSE)
}
