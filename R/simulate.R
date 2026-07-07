# ============================================================
# simulate.R
# simulate_model() -- run a FULLY-SPECIFIED compartmental sheet forward with no
# fitting and no observed data. Every state/parameter is fixed (*name=value), so
# there is nothing to estimate: the model is simply solved over the time grid and
# any dataDummy formulas are evaluated on the resulting trajectory.
#
#   sim <- simulate_model(modelParams, data_dummy = dummy)
#   sim$sir_out      # the solved trajectory (one column per compartment + date)
#   sim$evaluation   # each dummy/overlay formula evaluated over the grid
#   plot_simulation(sim)$grid
#
# This is the deterministic counterfactual mechanism (build -> recover the fixed
# point -> solve/evaluate) promoted to a first-class, data-free entry point.
# ============================================================

#' Simulate a fully-specified compartmental model (no fitting)
#'
#' Solve a compartmental ODE model whose states and parameters are ALL fixed
#' (`*name=value`), with no observed data and no estimation. The model is built,
#' its fixed point recovered, solved over the time grid, and any `data_dummy`
#' (and optional `dataCombined`) formulas are evaluated on the trajectory.
#'
#' Use this to run a scenario forward from known values -- e.g. a hand-specified
#' SIR/SEIR model -- without supplying `dataCombined`. If any quantity is still
#' fittable (a `[lo,hi]` box or a distributional prior), the sheet is not fully
#' specified and an error is raised naming the offending quantities; fit those
#' with [fitCompartmentalModel()] instead.
#'
#' @param modelParams Model parameter sheet (data frame). Every state and
#'   parameter must be fixed.
#' @param data_dummy Optional dummy-data data frame (a `Formula` column of
#'   expressions to evaluate on the trajectory; no likelihood, never fit).
#' @param dataCombined Optional observed data to also evaluate as overlay series
#'   (never fit). `NULL` (the default) simulates with dummy series only.
#' @param solver Solver settings, see [solver_control()]. Use
#'   `solver_control(backend = "r")` for a Julia-free solve.
#' @return An object of class `"compartmentalSim"`: a list with `initial_state`,
#'   `parms`, `sir_out` (solved trajectory), `evaluation` (formulas over the
#'   grid), `time_grid`, `model`, `solver`, `data_dummy` and `dataCombined`.
#' @seealso [plot_simulation()], [build_compartmental_model()],
#'   [fitCompartmentalModel()].
#' @examples
#' sim_dir <- system.file("extdata", "SIR_sim", package = "compfit")
#' mp <- read_data_file(file.path(sim_dir, "modelParams.csv"))
#' dd <- read_data_file(file.path(sim_dir, "dataDummy.csv"))
#' sim <- simulate_model(mp, data_dummy = dd, solver = solver_control(backend = "r"))
#' head(sim$sir_out)
#' sim$evaluation
#' @export
simulate_model <- function(modelParams,
                           data_dummy   = NULL,
                           dataCombined = NULL,
                           solver       = solver_control()) {
  # Build the machinery (validates the sheet, registers the ODE). With no
  # dataCombined this builds in simulation mode (loss = NULL).
  mod <- build_compartmental_model(modelParams, dataCombined, solver = solver)
  sap <- mod$sap

  # Require a FULLY-FIXED sheet: nothing left to estimate. Checked before the
  # zero-length recover call so the user gets an actionable message rather than
  # the order-contract guard tripping deep inside .recover_solution().
  fittable <- c(names(sap$states_fitted), names(sap$params_fitted))
  if (length(fittable))
    stop(sprintf(
      paste0("simulate_model() needs a fully-specified sheet, but these ",
             "quantities are still fittable: %s.\n  Fix each with *name=value, ",
             "or use fitCompartmentalModel() to estimate them."),
      paste(fittable, collapse = ", ")), call. = FALSE)

  # Recover every state/parameter from the fixed sheet (zero fitted quantities),
  # then solve and evaluate the dummy/overlay formulas on the trajectory.
  pe  <- .recover_solution(setNames(numeric(0), character(0)), sap, mod$time_grid)
  sae <- solve_and_evaluate(mod, pe$initial_state, pe$parms, data_dummy)

  structure(list(
    initial_state = pe$initial_state,
    parms         = pe$parms,
    sir_out       = sae$sir_out,
    evaluation    = sae$evaluation,
    time_grid     = mod$time_grid,
    model         = mod$model,
    solver        = solver,
    data_dummy    = data_dummy,
    dataCombined  = dataCombined
  ), class = "compartmentalSim")
}

#' @export
print.compartmentalSim <- function(x, ...) {
  tg <- x$time_grid
  cat("<compartmentalSim>\n")
  cat(sprintf("  compartments : %s\n", paste(names(x$initial_state), collapse = ", ")))
  cat(sprintf("  time span    : %d-%d (partition %d)\n",
              tg$startpoint, tg$endpoint, tg$partition))
  cat(sprintf("  backend      : %s\n",
              if (is.null(x$solver$backend)) "julia" else x$solver$backend))
  ps <- x$parms
  if (length(ps)) {
    shown <- paste(sprintf("%s=%g", names(ps), as.numeric(ps)), collapse = ", ")
    cat(sprintf("  parameters   : %s\n", shown))
  }
  series <- setdiff(names(x$evaluation), "date")
  if (length(series))
    cat(sprintf("  evaluated    : %s\n", paste(series, collapse = ", ")))
  cat(sprintf("  trajectory   : %d time points\n", nrow(x$sir_out)))
  invisible(x)
}

# ---- Plotting --------------------------------------------------------------

#' Plot a simulated model
#'
#' Plot the trajectories from a [simulate_model()] result: one panel per
#' compartment and/or per evaluated formula (dummy / overlay series), each a
#' line over the date grid. Mirrors [plot_fit()]'s return shape.
#'
#' @param sim A `"compartmentalSim"` object from [simulate_model()].
#' @param which Which panels to draw: `"states"` (compartment trajectories),
#'   `"series"` (evaluated dummy/overlay formulas), or `"both"` (default).
#' @param palette Colour palette: a scheme name (`"okabe"`/`"grey"`), a custom
#'   list like [cfit_palette], or `NULL` to use the `compfit.palette` option.
#' @param ncol Number of columns in the arranged grid.
#' @param base_size Base font size for the theme.
#' @return A list with `plots` (named list of ggplot objects) and `grid` (a
#'   patchwork arrangement, or `NULL` if patchwork is unavailable).
#' @seealso [simulate_model()], [plot_fit()].
#' @examples
#' sim_dir <- system.file("extdata", "SIR_sim", package = "compfit")
#' mp <- read_data_file(file.path(sim_dir, "modelParams.csv"))
#' dd <- read_data_file(file.path(sim_dir, "dataDummy.csv"))
#' sim <- simulate_model(mp, data_dummy = dd, solver = solver_control(backend = "r"))
#' \dontrun{ plot_simulation(sim)$grid }
#' @export
plot_simulation <- function(sim,
                            which     = c("both", "states", "series"),
                            palette   = NULL,
                            ncol      = 2,
                            base_size = 11) {
  if (!inherits(sim, "compartmentalSim"))
    stop("plot_simulation() expects a 'compartmentalSim' from simulate_model().")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_simulation() needs the 'ggplot2' package (install.packages(\"ggplot2\")).")
  which <- match.arg(which)
  pal   <- .cfit_resolve_palette(palette)

  # Assemble (series-name, source-data-frame, colour) panels.
  panels <- list()
  if (which %in% c("both", "states")) {
    for (nm in names(sim$initial_state))
      panels[[nm]] <- list(df = sim$sir_out, col = nm, colour = pal$model)
  }
  if (which %in% c("both", "series")) {
    # A dummy/overlay formula may coincide with a compartment name (e.g. "X1");
    # in "both" mode the compartment panel already covers it, so don't duplicate.
    for (nm in setdiff(names(sim$evaluation), "date"))
      if (!nm %in% names(panels))
        panels[[nm]] <- list(df = sim$evaluation, col = nm, colour = pal$dummy)
  }
  if (!length(panels))
    stop("Nothing to plot: no compartments or evaluated series selected.")

  brks_of <- function(dates) {
    yrs <- sort(unique(as.integer(format(dates, "%Y"))))
    as.Date(sprintf("%d-12-31", yrs[seq(1, length(yrs), length.out = min(6, length(yrs)))]))
  }
  plots <- lapply(names(panels), function(nm) {
    pv <- panels[[nm]]
    d  <- data.frame(date = pv$df$date, y = pv$df[[pv$col]])
    ggplot2::ggplot(d, ggplot2::aes(x = date, y = y)) +
      ggplot2::geom_line(colour = pv$colour, linewidth = 0.8) +
      ggplot2::scale_x_date(breaks = brks_of(d$date), date_labels = "%Y") +
      # Explicit y-axis number format so every panel is consistent (and matches
      # plot_fit): a thousands comma, never scientific. Without this the ggplot
      # default labeller is version-dependent and can add commas on some panels
      # but not others (e.g. "1,000" for one compartment, "5000" for another).
      ggplot2::scale_y_continuous(
        n.breaks = 5,
        labels = function(v) format(v, big.mark = ",", scientific = FALSE,
                                    trim = TRUE)) +
      ggplot2::labs(title = nm, x = NULL, y = NULL) +
      cfit_theme(base_size = base_size)
  })
  names(plots) <- names(panels)

  grid <- NULL
  if (requireNamespace("patchwork", quietly = TRUE))
    grid <- patchwork::wrap_plots(plots, ncol = ncol)

  list(plots = plots, grid = grid)
}
