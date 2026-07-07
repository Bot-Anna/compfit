lossFunction <- function(names_data_points,
                         parms_expression,
                         initial_states_expression,
                         states_params_expression,
                         penalty_expression,
                         time,
                         startpoint,
                         endpoint,
                         partition,
                         matrix_data_points,
                         weight_matrix,
                         average_matrix,
                         obs_mask,
                         cens_mask,
                         limit_mat,
                         lcens_mask,
                         llimit_mat,
                         interval_mask = NULL,
                         ilow_mat      = NULL,
                         iupp_mat      = NULL,
                         idev_lo_mat   = NULL,
                         idev_hi_mat   = NULL,
                         asym_mask     = NULL,
                         asym_val_mat  = NULL,
                         asym_dev_mat  = NULL,
                         asym_dir_mat  = NULL,
                         lower_guesses,
                         upper_guesses,
                         comp_names = NULL,
                         verbose = TRUE,
                         solver  = "AutoTsit5(Rosenbrock23())",
                         abstol  = 1e-8,
                         reltol  = 1e-8,
                         backend = "julia",
                         checkpoint_file = "checkpoint_best_solution.rds"
) {
  
  ## ---- Marker-based body insertion ----
  # Splices a parsed code block in place of a bare-symbol marker statement in
  # f's body. Parses ONLY code_str (never re-deparses the existing body), and
  # matches the marker as a language object rather than by text.
  # Requires each marker to appear as its own statement exactly once.
  insert_at_marker <- function(f, marker, code_str) {
    b <- as.list(body(f))                           # body is `{ ... }`; b[[1]] is `{`
    marker_sym <- as.symbol(marker)
    pos <- which(vapply(b, identical, logical(1), marker_sym))
    if (length(pos) != 1L)
      stop(sprintf("Marker '%s' not found exactly once.", marker))
    inserted <- as.list(parse(text = code_str))     # parsed statements, never re-deparsed
    b <- append(b[-pos], inserted, after = pos - 1) # drop marker, splice block in its place
    body(f) <- as.call(b)
    f
  }
  
  ## ---- We build the error expression for within the loss_function ----
  # Define integrate_for_loss inside lossFunction scope
  integrate_block_fn_str <- "
      integrate_for_loss <- function(f_vec, t_vec, part, start, end) {
      n   <- part * (end - start + 1)
      grp <- rep(seq_len(end - start + 1), each = part)
      as.numeric(tapply(seq_len(n), grp,
        function(idx) pracma::trapz(t_vec[idx], as.numeric(f_vec)[idx])))
      }
  "
  
  # Rewrite compartment references in a data formula to solution columns. A
  # compartment may be referenced by its X<i> position (identity case) or by its
  # declared name (S, I, R, ...); either way it maps to sir_out_data$<name>.
  # Longest names first so a short name can't clobber a longer one it prefixes.
  .sub_states <- function(s) {
    if (is.null(comp_names) || !length(comp_names))
      return(gsub("\\b(X\\d+)\\b", "sir_out_data$\\1", s))
    for (nm in comp_names[order(-nchar(comp_names))])
      s <- gsub(paste0("\\b\\Q", nm, "\\E\\b"), paste0("sir_out_data$", nm), s, perl = TRUE)
    s
  }

  make_integrate_block <- function(inner_str) {
    inner_str <- .sub_states(inner_str)
    paste0(
      "integrate_for_loss(",
      inner_str, ", ",
      "sir_out_data$time, partition, startpoint, endpoint)"
    )
  }
  
  error_expression <- paste0(
    integrate_block_fn_str,
    "\nannual_idx_loss <- seq(partition + 1, partition * (endpoint - startpoint + 1) + 1, by = partition)",
    "\nmatrix_sir <- cbind("
  )
  
  for (i in seq_along(names_data_points)) {
    names_expr <- gsub("`", "", names_data_points[i])
    
    if (grepl("^annual\\(", names_expr)) {
      inner_str   <- sub("^annual\\((.*)\\)$", "\\1", names_expr)
      expr_string <- make_integrate_block(inner_str)
    } else if (grepl("^cumulative\\(", names_expr)) {
      inner_str   <- sub("^cumulative\\((.*)\\)$", "\\1", names_expr)
      inner_str_r <- .sub_states(inner_str)
      # Cumulative integral at each annual snapshot
      expr_string <- paste0(
        "(function() {",
        " flux <- as.numeric(", inner_str_r, ");",
        " dt <- diff(sir_out_data$time);",
        " cumflux <- c(0, cumsum(0.5 * (flux[-length(flux)] + flux[-1]) * dt));",
        " cumflux[annual_idx_loss]",
        "})()"
      )
    } else {
      expr_string <- .sub_states(names_expr)
      expr_string <- paste0("(", expr_string, ")[annual_idx_loss]")
    }
    
    error_expression <- paste0(error_expression,
                               "\nV", i, " = ", expr_string, ",")
  }
  error_expression <- substr(error_expression, 1, nchar(error_expression) - 1)
  error_expression <- paste(c(error_expression, ")\n"), collapse = "")
  
  ## ---- Default the interval/asymmetric matrices (absent for older fits) ----
  # Coerce NULL -> a neutral matrix matching the data so the closure (which
  # captures these) can reference them unconditionally. dev defaults to 1 to
  # avoid a divide-by-zero when the asymmetric block is skipped.
  .z0 <- matrix(0, nrow(matrix_data_points), ncol(matrix_data_points))
  if (is.null(interval_mask)) interval_mask <- .z0
  if (is.null(ilow_mat))      ilow_mat      <- .z0
  if (is.null(iupp_mat))      iupp_mat      <- .z0
  if (is.null(idev_lo_mat))   idev_lo_mat   <- .z0
  if (is.null(idev_hi_mat))   idev_hi_mat   <- .z0
  if (is.null(asym_mask))     asym_mask     <- .z0
  if (is.null(asym_val_mat))  asym_val_mat  <- .z0
  if (is.null(asym_dir_mat))  asym_dir_mat  <- .z0
  if (is.null(asym_dev_mat))  asym_dev_mat  <- .z0 + 1

  ## ---- Set up of loss function ----
  ### Markers (.PARMS_PLACEHOLDER. etc.) are bare symbols, each on its own line,
  ### replaced below via insert_at_marker.
  loss_function <- function(states_and_params,
                            best_state,
                            show_error = FALSE){
    
    # The spliced blocks below look up parameters/states BY NAME in
    # states_and_params. The names come from lower_guesses (the fitted-quantity
    # order), so we attach them here rather than relying on the caller to pass a
    # named vector. This makes the loss callable with a bare numeric vector.
    states_and_params <- setNames(as.numeric(states_and_params),
                                  names(lower_guesses))
    states_and_params <- mapply(denormalise, states_and_params,
                                lower_guesses, upper_guesses)
    
    # Define the per-state/param _0 variables FIRST; parms and initial_states
    # below are built from them, so this block must precede both.
    .STATES_PARAMS_PLACEHOLDER.
    .PARMS_PLACEHOLDER.
    .INITIAL_STATES_PLACEHOLDER.
    
    ## ---- Julia part ----
    X <- unlist(initial_states)
    p <- unlist(parms)
    t <- c(as.numeric(min(time)), as.numeric(max(time)))
    .SOLVE_CALL.
    
    sol_matrix <- sol_result$matrix
    sir_out <- as.data.frame(t(sol_matrix))
    colnames(sir_out) <- names(initial_states)
    sir_out$time <- sol_result$t
    
    sir_out_data <- as.data.frame(sir_out)
    
    # Make each parameter available as a bare variable (e.g. `g_HA`) in addition
    # to the `parms` vector, so data formulas can be written as `g_HA*...`
    # instead of `parms["g_HA"]*...`. Both styles work (backward compatible).
    for (.nm in names(parms)) assign(.nm, parms[[.nm]])
    
    .ERROR_PLACEHOLDER.
    
    # Build the residual with censoring/missing handling:
    #  - observed (obs_mask==1): ordinary residual (matrix_sir - data)
    #  - censored (cens_mask==1): one-sided penalty max(0, model - L); no penalty
    #    while the model is at or below the limit, penalty if it exceeds it.
    #  - missing (neither): residual forced to 0 (contributes nothing).
    # Guard a failed / blown-up solve: non-finite model values would be zeroed
    # below (is.na(NaN) == TRUE), making a failure look like a PERFECT fit and
    # luring the optimiser into the failure region. Return a large finite penalty
    # instead, so the optimiser avoids it without seeing an Inf/NaN objective.
    if (!all(is.finite(sol_matrix)) || !all(is.finite(matrix_sir)))
      return(1e12)

    resid <- matrix_sir - matrix_data_points
    resid[is.na(resid)] <- 0                       # genuinely missing data (NA) -> 0
    if (any(cens_mask == 1)) {
      over <- pmax(0, matrix_sir - limit_mat)      # exceedance above upper bound
      over[is.na(over)] <- 0
      resid[cens_mask == 1] <- over[cens_mask == 1]
    }
    if (any(lcens_mask == 1)) {
      under <- pmax(0, llimit_mat - matrix_sir)    # shortfall below lower bound
      under[is.na(under)] <- 0
      resid[lcens_mask == 1] <- under[lcens_mask == 1]
    }
    if (any(interval_mask == 1)) {
      # interval [A,B]: two-sided hinge -- zero inside, linear outside. Soft
      # shoulders [A,B]~s divide each side by its scale (0 => weight 1 = hard).
      wlo <- ifelse(idev_lo_mat > 0, 1 / idev_lo_mat, 1)
      whi <- ifelse(idev_hi_mat > 0, 1 / idev_hi_mat, 1)
      outside <- pmax(0, ilow_mat - matrix_sir) * wlo + pmax(0, matrix_sir - iupp_mat) * whi
      outside[is.na(outside)] <- 0
      resid[interval_mask == 1] <- outside[interval_mask == 1]
    }
    if (any(asym_mask == 1)) {
      # asymmetric: hard hinge against the soft direction + soft (/dev) cost with it.
      hard_below <- pmax(0, asym_val_mat - matrix_sir)   # A - mu, clipped
      hard_above <- pmax(0, matrix_sir - asym_val_mat)   # mu - A, clipped
      up_pen   <- hard_below + hard_above / asym_dev_mat # soft upward   (dir = +1)
      down_pen <- hard_above + hard_below / asym_dev_mat # soft downward (dir = -1)
      asym_pen <- ifelse(asym_dir_mat > 0, up_pen, down_pen)
      asym_pen[is.na(asym_pen)] <- 0
      resid[asym_mask == 1] <- asym_pen[asym_mask == 1]
    }
    resid[obs_mask == 0 & cens_mask == 0 & lcens_mask == 0 &
          interval_mask == 0 & asym_mask == 0] <- 0  # missing
    
    error <- norm(resid %*% weight_matrix %*% average_matrix, type = "F")^2
    
    .PENALTY_PLACEHOLDER.
    
    if (show_error == TRUE) {
      if (error < best_state$error) {
        # Progress printing is OPTIONAL (verbose); best-state tracking +
        # checkpointing below happen regardless, so a quiet fit still recovers
        # the best solution.
        if (isTRUE(verbose)) print(error)
        best_state$error <- unname(error)
        best_state$par   <- states_and_params
        
        saveRDS(
          list(
            error = best_state$error,
            par   = best_state$par,
            time  = Sys.time()
          ),
          file = checkpoint_file
        )
      }
    }
    
    return(error)
  }
  
  # We splice the generated code blocks in at their markers.
  loss_function <- insert_at_marker(loss_function, ".PARMS_PLACEHOLDER.",          parms_expression)
  loss_function <- insert_at_marker(loss_function, ".INITIAL_STATES_PLACEHOLDER.", initial_states_expression)
  loss_function <- insert_at_marker(loss_function, ".STATES_PARAMS_PLACEHOLDER.",  states_params_expression)
  loss_function <- insert_at_marker(loss_function, ".ERROR_PLACEHOLDER.",          error_expression)
  loss_function <- insert_at_marker(loss_function, ".PENALTY_PLACEHOLDER.",        penalty_expression)
  
  # Solver settings are constants for the life of this loss closure; bake them
  # in. The backend decides which solver the loss calls: Julia (fast) or the
  # pure-R deSolve path (no Julia required).
  solve_call <- if (identical(backend, "r")) {
    sprintf(
      'sol_result <- solveWithR(X, t, p, time, abstol = %g, reltol = %g, method = "%s")',
      abstol, reltol, .desolve_method(solver)
    )
  } else {
    sprintf(
      'sol_result <- solveWithJulia(X, t, p, time, solver = "%s", abstol = %g, reltol = %g)',
      solver, abstol, reltol
    )
  }
  loss_function <- insert_at_marker(loss_function, ".SOLVE_CALL.", solve_call)
  
  return(loss_function)
}