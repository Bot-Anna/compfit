# ============================================================
# evaluate.R
# Solve the ODE at a given (initial_state, parms) and evaluate every data and
# dummy formula on the resulting trajectory. Self-sufficient: reads the time
# grid, solver settings, formulas and dummies from the `fit` object, so it needs
# no global variables. Used by main.R for plotting and by plot_fit().
#
#   ev <- solve_and_evaluate(fit, initial_state, parms)
#   ev$sir_out      # raw trajectory (states by time) + date column
#   ev$evaluation   # one column per data/dummy formula + date
#
# `data_dummy` is optional; pass NULL (default) to evaluate only data streams.
# ============================================================

#' Solve the ODE and evaluate formulas
#'
#' Solves the model at a given `(initial_state, parms)` and evaluates every data
#' (and optional dummy) formula on the resulting trajectory. Self-sufficient:
#' reads the time grid, solver settings and formulas from the `fit` object.
#'
#' @param fit A `"compartmentalFit"` object.
#' @param initial_state Named numeric vector of initial compartment values.
#' @param parms Named numeric vector of parameters.
#' @param data_dummy Optional dummy-data data frame; `NULL` evaluates only the
#'   data streams.
#' @return A list with `sir_out` (trajectory) and `evaluation` (formula columns).
#' @examples
#' \dontrun{
#' # fit from fitCompartmentalModel(); evaluate at the fitted point
#' p  <- get_point(fit)
#' ev <- solve_and_evaluate(fit, p$initial_state, p$parms)
#' head(ev$evaluation)
#' }
#' @export
solve_and_evaluate <- function(fit, initial_state, parms, data_dummy = NULL) {
  tg        <- fit$time_grid
  time      <- tg$time
  partition <- tg$partition          # passed explicitly to evaluate_formula()
  date      <- fit$model$date
  
  X <- unlist(initial_state)
  p <- unlist(parms)
  t <- c(as.numeric(min(time)), as.numeric(max(time)))
  
  sol <- if (identical(.ode_backend(fit$solver$backend), "r")) {
    solveWithR(X, t, p, time,
               abstol = fit$solver$abstol, reltol = fit$solver$reltol,
               method = .desolve_method(fit$solver$solver))
  } else {
    solveWithJulia(X, t, p, time,
                   solver = fit$solver$solver,
                   abstol = fit$solver$abstol,
                   reltol = fit$solver$reltol)
  }
  sir_out <- as.data.frame(t(sol$matrix))
  colnames(sir_out) <- names(initial_state)
  
  # Guard: a failed/aborted ODE solve returns a trajectory whose length does not
  # match the time grid. Catch it here with a legible message rather than letting
  # the date-column assignment fail cryptically ("replacement has N rows ...").
  if (nrow(sir_out) != length(date)) {
    stop(sprintf(
      "ODE solve returned %d time points but the grid has %d. The solve likely ",
      nrow(sir_out), length(date)),
      "failed/aborted at these parameter values (often a stiff/unstable region, ",
      "e.g. plotting at a non-converged posterior mean). Check the parameter ",
      "values and the solver warnings above.")
  }
  sir_out$date <- date
  
  streams <- fit$data$names_data_points
  ev <- data.frame(date = sir_out$date)
  for (col in streams) {
    ev[[col]] <- evaluate_formula(col, sir_out, parms, time, partition)
  }
  if (!is.null(data_dummy)) {
    for (col in data_dummy$Formula) {
      ev[[col]] <- evaluate_formula(col, sir_out, parms, time, partition)
    }
  }
  
  list(sir_out = sir_out, evaluation = ev)
}