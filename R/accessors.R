# ============================================================
# accessors.R
# A stable, discoverable accessor layer over a `compartmentalFit` object.
#
# The fit object's fields remain public (you can still write fit$loss etc.),
# but these accessors give a documented, validated API that won't break if the
# internal layout shifts, and they fail with a clear message when something
# isn't available (e.g. asking for posterior samples on an MLE fit).
#
#   get_loss(fit)                  the loss closure
#   get_model(fit)                 the model bundle (structure/expressions/etc.)
#   get_compartmental_function(fit)the R compartmental function
#   get_julia_code(fit)            the registered Julia ODE source
#   get_expressions(fit)           generated expression pieces
#   get_states_and_params(fit)     the statesAndParams() output (sap)
#   get_data(fit)                  the prepared-data bundle
#   get_time_grid(fit)             list(time, startpoint, endpoint, partition, cutoff)
#   get_bounds(fit)                list(lower, upper, init_norm)
#   get_solver(fit)                solver_control() used
#   get_point(fit)                 list(initial_state, parms)  [MLE-type fits]
#   get_samples(fit)               raw Bayes sampler return
# ============================================================

.is_fit <- function(fit) {
  if (!inherits(fit, "compartmentalFit"))
    stop("Expected a 'compartmentalFit' object (the result of fitCompartmentalModel()).")
  invisible(TRUE)
}

#' Accessors for a compartmental fit
#'
#' Extract components of a `"compartmentalFit"` object.
#'
#' @param fit A `"compartmentalFit"` object.
#' @return The requested component of the fit.
#' @name compfit-accessors
#' @examples
#' \dontrun{
#' # fit from fitCompartmentalModel()
#' get_point(fit)
#' get_julia_code(fit)
#' }
#' @export
get_loss        <- function(fit) { .is_fit(fit); fit$loss }
#' @rdname compfit-accessors
#' @export
get_model       <- function(fit) { .is_fit(fit); fit$model }
#' @rdname compfit-accessors
#' @export
get_data        <- function(fit) { .is_fit(fit); fit$data }
#' @rdname compfit-accessors
#' @export
get_time_grid   <- function(fit) { .is_fit(fit); fit$time_grid }
#' @rdname compfit-accessors
#' @export
get_bounds      <- function(fit) { .is_fit(fit); fit$bounds }
#' @rdname compfit-accessors
#' @export
get_solver      <- function(fit) { .is_fit(fit); fit$solver }

#' @rdname compfit-accessors
#' @export
get_compartmental_function <- function(fit) {
  .is_fit(fit); fit$model$compartmental_function
}
#' @rdname compfit-accessors
#' @export
get_julia_code  <- function(fit) { .is_fit(fit); fit$model$julia_code }
#' @rdname compfit-accessors
#' @export
get_expressions <- function(fit) { .is_fit(fit); fit$model$expressions }
#' @rdname compfit-accessors
#' @export
get_states_and_params <- function(fit) { .is_fit(fit); fit$model$sap }

#' @rdname compfit-accessors
#' @export
get_point <- function(fit) {
  .is_fit(fit)
  if (is.null(fit$point))
    stop("This fit has no point estimate (it is a Bayesian fit, or the fit failed). ",
         "Use get_samples(fit) / posterior_*(fit) for a Bayesian fit.")
  fit$point
}

#' @rdname compfit-accessors
#' @export
get_samples <- function(fit) {
  .is_fit(fit)
  if (is.null(fit$samples))
    stop("This fit has no posterior samples (it is not a method='bayes' fit, or it failed). ",
         "Use get_point(fit) for an MLE-type fit.")
  fit$samples
}

# ---- summary method: what's in this fit and how to reach it ----------------
#' Summary of a compartmental fit
#'
#' @param object A `"compartmentalFit"` object.
#' @param ... Unused.
#' @return `object`, invisibly (called for the printed summary).
#' @examples
#' \dontrun{
#' summary(fit)   # fit from fitCompartmentalModel()
#' }
#' @exportS3Method summary compartmentalFit
summary.compartmentalFit <- function(object, ...) {
  x <- object
  cat(sprintf("compartmentalFit  |  method = %s  |  success = %s\n",
              x$method, x$success))
  if (!isTRUE(x$success) && !is.null(x$error_msg))
    cat("  error:", x$error_msg, "\n")
  
  tg <- x$time_grid
  cat(sprintf("  horizon: %d..%d  (%d years, partition = %d)\n",
              tg$startpoint, tg$endpoint, tg$endpoint - tg$startpoint + 1, tg$partition))
  cat(sprintf("  compartments: %d   data streams: %d\n",
              x$model$structure$number_of_comps, length(x$data$names_data_points)))
  
  if (x$method == "bayes" && !is.null(x$samples)) {
    d <- x$samples$draws
    cat(sprintf("  posterior: %d draws x %d sampled quantities\n", nrow(d), ncol(d)))
    cat("  access: get_samples(fit), posterior_report(fit), posterior_summary(fit)\n")
  } else if (!is.null(x$point)) {
    cat("  point estimate: get_point(fit)$parms / $initial_state\n")
  }
  
  cat("  objects: get_loss(), get_model(), get_compartmental_function(),\n",
      "          get_julia_code(), get_expressions(), get_states_and_params(),\n",
      "          get_data(), get_time_grid(), get_bounds(), get_solver()\n", sep = "")
  invisible(x)
}