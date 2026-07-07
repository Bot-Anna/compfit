hypercubeSampling <- function(n, d, loss_function,
                              not_sobol_states_normalised,
                              not_sobol_params_normalised,
                              lower_guesses_states_without,
                              upper_guesses_states_without,
                              lower_guesses_params_without,
                              upper_guesses_params_without,
                              lower_guesses,
                              upper_guesses,
                              best_state,
                              progress = TRUE) {
  number_of_samples     <- n
  dimension_of_hypercube <- d
  sobol_points          <- sobol(n = number_of_samples, d = dimension_of_hypercube)
  center                <- rep(0.5, dimension_of_hypercube)
  sobol_points_centred  <- (sobol_points + matrix(center,
                                                  nrow = number_of_samples,
                                                  ncol = dimension_of_hypercube,
                                                  byrow = TRUE)) %% 1
  
  best_error  <- Inf
  best_sol    <- vector()

  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = number_of_samples, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in 1:number_of_samples) {
    if (isTRUE(progress)) utils::setTxtProgressBar(pb, i)
    # seq_len (not 1:length): when there are no "without" states, 1:length(0)
    # would be c(1, 0) and grab column 1 instead of nothing. seq_len(0) is empty.
    states_without <- sobol_points_centred[i, seq_len(length(lower_guesses_states_without))]
    names(states_without) <- names(lower_guesses_states_without)
    
    params_without <- sobol_points_centred[i, (length(lower_guesses_states_without) + 1):dimension_of_hypercube]
    names(params_without) <- names(lower_guesses_params_without)
    
    initial_guess_sobol <- c(not_sobol_states_normalised,
                             states_without,
                             not_sobol_params_normalised,
                             params_without)
    initial_guess_sobol <- initial_guess_sobol[!is.na(names(initial_guess_sobol))]
    
    current_error <- tryCatch(
      loss_function(initial_guess_sobol, best_state, show_error = FALSE),
      error = function(e) Inf
    )
    
    if (is.finite(current_error) && current_error < best_error) {
      if (!isTRUE(progress))   # the progress bar already shows activity
        cat("Sample", i, "- new best:", unname(current_error), "\n")
      best_error <- current_error
      best_sol <- mapply(denormalise,
                         initial_guess_sobol,
                         lower_guesses,
                         upper_guesses)
      best_state$error <- unname(current_error)
      # Store the NATURAL-scale solution (not the normalised initial_guess_sobol):
      # .recover_solution()/the loss store the denormalised par here, so the
      # downstream point recovery in fitCompartmentalModel() works for hypercube
      # fits exactly as it does for optim/deoptim.
      best_state$par   <- best_sol
    }
  }
  
  return(list(best_sol = best_sol, best_state = best_state))
}