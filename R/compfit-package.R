#' compfit: Fit, Analyse and Compare Compartmental ODE Models (R + Julia)
#'
#' Build compartmental (SIR-type) ODE models from a spreadsheet specification,
#' fit them by maximum likelihood (L-BFGS-B / DEoptim / hypercube) or Bayesian
#' sampling, and analyse the results: posterior summaries, prior-vs-posterior
#' and predictive plots, global (Sobol) and local (derivative) sensitivity,
#' practical identifiability, counterfactual scenarios, and self-contained code
#' export. ODE solves and Bayesian sampling run in Julia via the JuliaCall
#' bridge; call [setup_julia()] once per session before fitting.
#'
#' @keywords internal
#' @import stats
#' @import utils
#' @import ggplot2
#' @importFrom grDevices colorRampPalette dev.list dev.off
#' @importFrom pracma trapz
#' @importFrom qrng sobol
#' @importFrom DEoptim DEoptim DEoptim.control
"_PACKAGE"

# Silence R CMD check "no visible binding" NOTEs for non-standard evaluation:
# data-frame columns referenced bare in ggplot2 aes()/dplyr, placeholder symbols
# spliced into the loss closure's assembled expression, and tic/toc state.
utils::globalVariables(c(
  # loss-closure template placeholders (built via substitute/bquote)
  ".ERROR_PLACEHOLDER.", ".INITIAL_STATES_PLACEHOLDER.", ".PARMS_PLACEHOLDER.",
  ".PENALTY_PLACEHOLDER.", ".SOLVE_CALL.", ".STATES_PARAMS_PLACEHOLDER.",
  "initial_states", "parms", "sol_result", "matrix_sir", "vector_of_functions",
  # tic/toc timing state
  ".start_time",
  # NSE data-frame columns used in plot/report aes() and aggregation
  "absS", "draw", "e", "hi", "limit", "lo", "lower", "med", "output", "p1", "p2",
  "parameter", "part", "sensitivity", "total", "upper", "value", "y",
  # interval / asymmetric data-cell plot markers
  "low", "high", "dev", "dir"
))
