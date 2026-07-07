# ============================================================
# io.R  –  persist and restore a compartmentalFit object
#
#   save_fit(fit, path)   write to an RDS file
#   load_fit(path)        read back a fully functional fit object
#   reload_fit(path, register = TRUE)  load_fit() when register = TRUE; a
#                         results-only read (no Julia, no loss) when FALSE
#
# Design notes:
#   best_state  – environment replaced by a plain list (two scalars);
#                 reconstructed as a new env on load.
#   loss        – dropped entirely; all inputs are already stored in
#                 other fit fields, so it is rebuilt by load_fit().
#   model$compartmental_function – pure R closure, serialised as-is.
#
# load_fit() also re-registers the Julia ODE function (julia_code) in
# the running Julia session, which is required before loss can be called.
# ============================================================

#' Save a fit to disk
#'
#' Serialises a `"compartmentalFit"` to an `.rds` file (the non-portable loss
#' closure is dropped; [load_fit()] rebuilds it).
#'
#' @param fit A `"compartmentalFit"` object.
#' @param path Output path (`.rds` appended if missing).
#' @return The output path, invisibly.
#' @examples
#' \dontrun{
#' save_fit(fit, file.path(tempdir(), "fit.rds"))   # fit from fitCompartmentalModel()
#' }
#' @export
save_fit <- function(fit, path) {
  .is_fit(fit)

  if (!grepl("\\.rds$", path, ignore.case = TRUE))
    path <- paste0(path, ".rds")

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  s            <- fit
  s$best_state <- list(error = fit$best_state$error,
                       par   = fit$best_state$par)
  s$loss       <- NULL

  saveRDS(s, file = path)
  invisible(path)
}

#' Load a fit from disk
#'
#' Reads a `"compartmentalFit"` saved by [save_fit()] and rebuilds the loss
#' closure, re-registering the Julia ODE function in the running session.
#'
#' @param path Path to the `.rds` file (`.rds` appended if missing).
#' @return A `"compartmentalFit"` object.
#' @examples
#' \dontrun{
#' setup_julia()
#' fit <- load_fit(file.path(tempdir(), "fit.rds"))
#' }
#' @export
load_fit <- function(path) {
  if (!grepl("\\.rds$", path, ignore.case = TRUE))
    path <- paste0(path, ".rds")

  s <- readRDS(path)

  # Restore best_state as a mutable environment
  bs        <- new.env(parent = emptyenv())
  bs$error  <- s$best_state$error
  bs$par    <- s$best_state$par
  s$best_state <- bs

  # Re-register the ODE function for this fit's backend so the rebuilt loss can
  # solve (Julia function, or the R closure for backend = "r").
  if (identical(.ode_backend(s$solver$backend), "r")) {
    registerRODEFunction(s$model$compartmental_function)
  } else {
    registerJuliaODEFunction(s$model$julia_code)
  }

  # Rebuild the loss closure from the stored fit fields
  s$loss <- lossFunction(
    names_data_points         = s$data$names_data_points,
    parms_expression          = s$model$expressions$parms,
    initial_states_expression = s$model$expressions$initial_states,
    states_params_expression  = s$model$expressions$states_params,
    penalty_expression        = s$model$expressions$penalty,
    time                      = s$time_grid$time,
    startpoint                = s$time_grid$startpoint,
    endpoint                  = s$time_grid$endpoint,
    partition                 = s$time_grid$partition,
    matrix_data_points        = s$data$matrix_data_points,
    weight_matrix             = s$data$weight_matrix,
    average_matrix            = s$data$average_matrix,
    obs_mask                  = s$data$obs_mask,
    cens_mask                 = s$data$cens_mask,
    limit_mat                 = s$data$limit_mat,
    lcens_mask                = s$data$lcens_mask,
    llimit_mat                = s$data$llimit_mat,
    interval_mask             = s$data$interval_mask,   # NULL for pre-feature fits
    ilow_mat                  = s$data$ilow_mat,
    iupp_mat                  = s$data$iupp_mat,
    asym_mask                 = s$data$asym_mask,
    asym_val_mat              = s$data$asym_val_mat,
    asym_dev_mat              = s$data$asym_dev_mat,
    asym_dir_mat              = s$data$asym_dir_mat,
    lower_guesses             = s$bounds$lower,
    upper_guesses             = s$bounds$upper,
    comp_names                = s$sap$comp_names,   # NULL for pre-feature fits
    solver                    = s$solver$solver,
    abstol                    = s$solver$abstol,
    reltol                    = s$solver$reltol,
    backend                   = .ode_backend(if (is.null(s$solver$backend)) "julia" else s$solver$backend),
    checkpoint_file           = s$meta$checkpoint_file
  )

  class(s) <- "compartmentalFit"
  s
}

# reload_fit(): the entry point used by the counterfactual pipeline and by
# interactive reloads. With register = TRUE it returns a fully re-executable fit
# (re-registers the Julia ODE and rebuilds the loss via load_fit(), so the loss /
# solve_and_evaluate() / plot_fit() / Sobol all work). With register = FALSE it
# does a results-only read -- restores the object and best_state but does NOT
# touch Julia or rebuild the loss -- for inspecting draws / posterior_* without a
# Julia session.
#' Reload a fit (optionally re-registering Julia)
#'
#' Entry point used by the counterfactual pipeline and interactive reloads. With
#' `register = TRUE` returns a fully re-executable fit (re-registers the Julia
#' ODE and rebuilds the loss via [load_fit()]); with `register = FALSE` does a
#' results-only read (object + draws, no Julia session).
#'
#' @param path Path to the `.rds` file (`.rds` appended if missing).
#' @param register If `TRUE`, re-register Julia and rebuild the loss.
#' @return A `"compartmentalFit"` object.
#' @examples
#' \dontrun{
#' fit  <- reload_fit(file.path(tempdir(), "fit.rds"))                  # with Julia
#' draws <- reload_fit(file.path(tempdir(), "fit.rds"), register = FALSE) # results only
#' }
#' @export
reload_fit <- function(path, register = TRUE) {
  if (!grepl("\\.rds$", path, ignore.case = TRUE))
    path <- paste0(path, ".rds")

  if (isTRUE(register)) {
    if (!exists("registerJuliaODEFunction"))
      stop("reload_fit(register = TRUE) needs the package loaded first ",
           "(source setup.R) so the Julia bridge is available. Use ",
           "register = FALSE for results-only inspection.")
    return(load_fit(path))           # full restore: Julia + rebuilt loss
  }

  # Results-only: read the object, restore best_state as a mutable env, no Julia.
  s <- readRDS(path)
  bs        <- new.env(parent = emptyenv())
  bs$error  <- s$best_state$error
  bs$par    <- s$best_state$par
  s$best_state <- bs
  s$loss <- NULL                      # not rebuilt; needs the registered ODE
  class(s) <- "compartmentalFit"
  s
}
