generateExpressions <- function(number_of_comps,
                                states_fitted, states_fixed, states_functions,
                                params_fitted, params_fixed, params_functions,
                                conditions,
                                comp_names = paste0("X", seq_len(number_of_comps))) {
  # For setting up the loss function and the initial_states expression
  
  # First the states, then the params
  # We proceed as follows:
  # First we create the expressions for the fixed ones, followed by the 
  # initialisation of the model fits, then the ones that are functions of 
  # the others
  
  ## ---- Subroutine for SIR model function ----
  # Builds the string of the return call of the SIR model function
  return_expression <- "\nreturn(c("
  
  # For loops to set up the main expressions
  states_params_expression <- ""
  initial_states_expression <- "initial_states <- c("
  for (i in 1:number_of_comps) {
    # Compartment i is named comp_names[i] (X<i> in the identity case, or
    # S/I/R/... under a named registry). Its value slot is <name>_0, matching the
    # states_and_params[['<name>']] assignments built above. The names here become
    # the solution's column names, so data formulas can reference S/I/R directly.
    nm <- comp_names[i]
    initial_states_expression <- paste(c(initial_states_expression,
                                         nm,
                                         "=",
                                         nm,
                                         "_0",
                                         ","),
                                       collapse = "")
  }
  
  parms_expression <- "parms <- c("
  end_of_expression <- ")"
  sir_expression <- ""
  if (length(params_fitted) != 0) {
    for (i in 1:length(params_fitted)) {
      states_params_expression <- paste(c(states_params_expression, 
                                          paste(c("\n", 
                                                  unname(params_fitted[i]),
                                                  "_0",
                                                  "<-", 
                                                  "states_and_params[['", 
                                                  names(params_fitted)[i], 
                                                  "']]"), 
                                                collapse="")), 
                                        collapse="")
      parms_expression <- paste(c(parms_expression,
                                  names(params_fitted)[i],
                                  "=",
                                  names(params_fitted)[i],
                                  "_0",
                                  ","),
                                collapse = "")
      sir_expression <- paste(c(sir_expression, "\n",
                                names(params_fitted)[i],
                                "=p[",
                                i,
                                "]"),
                              collapse = "")
    }
  }
  if (length(states_fitted) != 0) {
    for (i in (length(params_fitted)+1):
      (length(states_fitted)+length(params_fitted))) {
      states_params_expression <- paste(c(states_params_expression, 
                                          paste(c("\n", 
                                                  unname(states_fitted[(i-length(params_fitted))]),
                                                  "_0",
                                                  "<-", 
                                                  "states_and_params[['", 
                                                  names(states_fitted)[(i-length(params_fitted))], 
                                                  "']]"), 
                                                collapse="")), 
                                        collapse="")
    }
  }
  if (length(params_fixed) != 0) {
    for (i in 1:length(params_fixed)) {
      states_params_expression <- paste(c(states_params_expression, 
                                          paste(c("\n", 
                                                  names(params_fixed)[i],
                                                  "_0",
                                                  "<-", 
                                                  unname(params_fixed[i])), 
                                                collapse="")), 
                                        collapse="")
      parms_expression <- paste(c(parms_expression,
                                  names(params_fixed)[i],
                                  "=",
                                  names(params_fixed)[i],
                                  "_0",
                                  ","),
                                collapse = "")
      sir_expression <- paste(c(sir_expression, "\n",
                                names(params_fixed)[i],
                                "=p[",
                                (i+length(params_fitted)),
                                "]"),
                              collapse = "")
    }
  }
  if (length(states_fixed) != 0) {
    for (i in 1:length(states_fixed)) {
      states_params_expression <- paste(c(states_params_expression, 
                                          paste(c("\n", 
                                                  names(states_fixed)[i],
                                                  "_0", 
                                                  "<-", 
                                                  unname(states_fixed[i])), 
                                                collapse = "")), 
                                        collapse = "")
    }
  }
  if (length(states_functions) != 0) {
    for (i in 1:length(states_functions)) {
      states_params_expression <- paste(c(states_params_expression, 
                                          paste(c("\n", 
                                                  names(states_functions)[i],
                                                  "_0",
                                                  "<-", 
                                                  unname(states_functions[i])), 
                                                collapse="")), 
                                        collapse="")
    }
  }
  if (length(params_functions) != 0) {
    for (i in 1:length(params_functions)) {
      states_params_expression <- paste(c(states_params_expression, 
                                          paste(c("\n", 
                                                  names(params_functions)[i],
                                                  "_0",
                                                  "<-", 
                                                  unname(params_functions[i])), 
                                                collapse="")), 
                                        collapse="")
      parms_expression <- paste(c(parms_expression,
                                  names(params_functions)[i],
                                  "=",
                                  0,
                                  ","),
                                collapse = "")
      sir_expression <- paste(c(sir_expression, "\n",
                                names(params_functions)[i], 
                                "=p[", 
                                (i+length(params_fitted)+length(params_fixed)),
                                "]"), 
                              collapse = "")
    }
  }
  substr(initial_states_expression, 
         nchar(initial_states_expression), 
         nchar(initial_states_expression)) <- end_of_expression
  substr(parms_expression, 
         nchar(parms_expression), 
         nchar(parms_expression)) <- end_of_expression
  
  for (i in 1:number_of_comps) {
    return_expression <- paste(c(return_expression, "dX", i, ","), collapse = "")
  }
  
  return_expression <- substr(return_expression, 1, nchar(return_expression)-1)
  return_expression <- paste(c(return_expression, "))"), collapse="")
  
  # We set up the penalty term obtained from the conditions set

  build_penalty_expression <- function(constraints,
                                       parms_names,
                                       states_names,
                                       penalty_scale = 1e8,
                                       approx_scale = 1e4) {
    if (length(constraints) == 0) return("")
    
    r_reserved <- c(
      "exp", "log", "log2", "log10", "sqrt", "abs", "sign",
      "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
      "sinh", "cosh", "tanh",
      "floor", "ceiling", "round", "trunc",
      "min", "max", "sum", "prod", "cumsum", "cumprod",
      "pmin", "pmax", "diff", "range",
      "Inf", "NaN", "NA", "TRUE", "FALSE", "pi", "e",
      "ifelse", "if_else", "is.na", "is.nan", "is.infinite"
    )
    
    classify_tokens <- function(expr) {
      # Protect numeric literals
      numerics <- regmatches(
        expr,
        gregexpr("(?<![A-Za-z_0-9])\\d+\\.?\\d*([eE][+-]?\\d+)?(?![A-Za-z_0-9])", expr, perl = TRUE)
      )[[1]]
      placeholders <- paste0("__NUM", seq_along(numerics), "__")
      for (i in seq_along(numerics)) {
        expr <- gsub(
          paste0("(?<![A-Za-z_0-9])", 
                 gsub("\\.", "\\\\.", numerics[i]),  # escape decimal point if present
                 "(?![A-Za-z_0-9])"),
          placeholders[i],
          expr,
          perl = TRUE
        )
      }
      
      tokens <- regmatches(
        expr,
        gregexpr("\\b[A-Za-z_][A-Za-z0-9_.]*\\b", expr)
      )[[1]]
      tokens <- setdiff(unique(tokens), c(r_reserved,
                                          grep("^__NUM\\d+__$", tokens, value = TRUE)))      
      # Classify each token
      from_parms  <- intersect(tokens, parms_names)
      from_states <- intersect(tokens, states_names)
      unknown     <- setdiff(tokens, c(parms_names, states_names))
      
      if (length(unknown) > 0) {
        warning(paste0(
          "Tokens not found in parms or states: ",
          paste(unknown, collapse = ", "),
          " -- left as-is"
        ))
      }
      
      list(
        tokens      = tokens,
        from_parms  = from_parms,
        from_states = from_states,
        unknown     = unknown,
        numerics    = numerics,
        placeholders = placeholders,
        expr_protected = expr
      )
    }
    
    substitute_vars <- function(expr) {
      info <- classify_tokens(expr)
      expr <- info$expr_protected
      
      # Substitute parms tokens
      for (tok in info$from_parms) {
        expr <- gsub(
          paste0("\\b", tok, "\\b"),
          paste0("parms['", tok, "']"),
          expr
        )
      }
      
      # Substitute states tokens
      for (tok in info$from_states) {
        expr <- gsub(
          paste0("\\b", tok, "\\b"),
          paste0("initial_states['", tok, "']"),
          expr
        )
      }
      
      # Restore numeric literals
      for (i in seq_along(info$numerics)) {
        expr <- sub(info$placeholders[i], info$numerics[i], expr, fixed = TRUE)
      }
      
      expr
    }
    
    substitute_vars_trajectory <- function(expr) {
      info <- classify_tokens(expr)
      expr <- info$expr_protected
      
      for (tok in info$from_parms) {
        expr <- gsub(
          paste0("\\b", tok, "\\b"),
          paste0("parms['", tok, "']"),
          expr
        )
      }
      
      for (tok in info$from_states) {
        expr <- gsub(
          paste0("\\b", tok, "\\b"),
          paste0("sir_out_data[['", tok, "']]"),
          expr
        )
      }
      
      for (i in seq_along(info$numerics)) {
        expr <- sub(info$placeholders[i], info$numerics[i], expr, fixed = TRUE)
      }
      
      expr
    }
    
    expressions <- sapply(constraints, function(con) {
      
      is_diff_trajectory <- grepl("^diff_trajectory:", con)
      con <- trimws(sub("^diff_trajectory:", "", con))
      is_trajectory <- grepl("^trajectory:", con) | is_diff_trajectory
      con <- trimws(sub("^trajectory:", "", con))
      is_endpoint <- grepl("^endpoint:", con)
      con <- trimws(sub("^endpoint:", "", con))
      
      if (grepl("~=", con)) {
        op <- "~="
      } else if (grepl("<=", con)) {
        op <- "<="
      } else if (grepl(">=", con)) {
        op <- ">="
      } else if (grepl("<", con)) {
        op <- "<"
      } else if (grepl(">", con)) {
        op <- ">"
      } else {
        stop(paste("No recognised operator in constraint:", con))
      }
      
      parts <- trimws(strsplit(con, op, fixed = TRUE)[[1]])
      lhs_r <- substitute_vars(parts[1])
      rhs_r <- substitute_vars(parts[2])
      
      diff_lr <- paste0("(", lhs_r, ") - (", rhs_r, ")")
      diff_rl <- paste0("(", rhs_r, ") - (", lhs_r, ")")
      
      switch(op,
             "~=" = {
               if (is_diff_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("diff(", lhs_s, ") - (", rhs_s, ")")
                 paste0("error <- error + ", approx_scale,
                        " * sum(((", diff_s, ") / (mean(abs(diff(", lhs_s, "))) + 1e-8))^2)")
               } else if (is_endpoint) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 paste0("error <- error + ", approx_scale,
                        " * ((tail(", lhs_s, ", 1) - tail(", rhs_s, ", 1))",
                        " / (abs(tail(", lhs_s, ", 1)) + 1e-8))^2")
               } else if (is_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("(", lhs_s, ") - (", rhs_s, ")")
                 paste0("error <- error + ", approx_scale,
                        " * sum(((", diff_s, ") / (mean(abs(", lhs_s, ")) + 1e-8))^2)")
               } else {
                 paste0("error <- error + ", approx_scale, " * (", diff_lr, ")^2")
               }
             },
             "<=" = {
               if (is_diff_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("diff(", lhs_s, ") - (", rhs_s, ")")
                 paste0("error <- error + ", penalty_scale,
                        " * sum(pmax(0, (", diff_s, ") / (mean(abs(diff(", lhs_s, "))) + 1e-8))^2)")
               } else if (is_endpoint) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 paste0("error <- error + ", penalty_scale,
                        " * pmax(0, (tail(", lhs_s, ", 1) - tail(", rhs_s, ", 1))",
                        " / (abs(tail(", lhs_s, ", 1)) + 1e-8))^2")
               } else if (is_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("(", lhs_s, ") - (", rhs_s, ")")
                 paste0("error <- error + ", penalty_scale,
                        " * sum(pmax(0, (", diff_s, ") / (mean(abs(", lhs_s, ")) + 1e-8))^2)")
               } else {
                 paste0("error <- error + ", penalty_scale, " * max(0, ", diff_lr, ")^2")
               }
             },
             ">=" = {
               if (is_diff_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("diff(", lhs_s, ") - (", rhs_s, ")")
                 paste0("error <- error + ", penalty_scale,
                        " * sum(pmax(0, -(", diff_s, ") / (mean(abs(diff(", lhs_s, "))) + 1e-8))^2)")
               } else if (is_endpoint) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 paste0("error <- error + ", penalty_scale,
                        " * pmax(0, (tail(", rhs_s, ", 1) - tail(", lhs_s, ", 1))",
                        " / (abs(tail(", lhs_s, ", 1)) + 1e-8))^2")
               } else if (is_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("(", rhs_s, ") - (", lhs_s, ")")
                 paste0("error <- error + ", penalty_scale,
                        " * sum(pmax(0, (", diff_s, ") / (mean(abs(", lhs_s, ")) + 1e-8))^2)")
               } else {
                 paste0("error <- error + ", penalty_scale, " * max(0, ", diff_rl, ")^2")
               }
             },
             "<" = {
               if (is_diff_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("diff(", lhs_s, ") - (", rhs_s, ")")
                 paste0("error <- error + ", penalty_scale,
                        " * sum(pmax(0, (", diff_s, ") / (mean(abs(diff(", lhs_s, "))) + 1e-8) + 1e-6)^2)")
               } else if (is_endpoint) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 paste0("error <- error + ", penalty_scale,
                        " * pmax(0, (tail(", lhs_s, ", 1) - tail(", rhs_s, ", 1))",
                        " / (abs(tail(", lhs_s, ", 1)) + 1e-8) + 1e-6)^2")
               } else if (is_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("(", lhs_s, ") - (", rhs_s, ")")
                 paste0("error <- error + ", penalty_scale,
                        " * sum(pmax(0, (", diff_s, ") / (mean(abs(", lhs_s, ")) + 1e-8) + 1e-6)^2)")
               } else {
                 paste0("error <- error + ", penalty_scale, " * max(0, ", diff_lr, " + 1e-6)^2")
               }
             },
             ">" = {
               if (is_diff_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("diff(", lhs_s, ") - (", rhs_s, ")")
                 paste0("error <- error + ", penalty_scale,
                        " * sum(pmax(0, -(", diff_s, ") / (mean(abs(diff(", lhs_s, "))) + 1e-8) + 1e-6)^2)")
               } else if (is_endpoint) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 paste0("error <- error + ", penalty_scale,
                        " * pmax(0, (tail(", rhs_s, ", 1) - tail(", lhs_s, ", 1))",
                        " / (abs(tail(", lhs_s, ", 1)) + 1e-8) + 1e-6)^2")
               } else if (is_trajectory) {
                 lhs_s <- substitute_vars_trajectory(parts[1])
                 rhs_s <- substitute_vars_trajectory(parts[2])
                 diff_s <- paste0("(", rhs_s, ") - (", lhs_s, ")")
                 paste0("error <- error + ", penalty_scale,
                        " * sum(pmax(0, (", diff_s, ") / (mean(abs(", lhs_s, ")) + 1e-8) + 1e-6)^2)")
               } else {
                 paste0("error <- error + ", penalty_scale, " * max(0, ", diff_rl, " + 1e-6)^2")
               }
             }
      )
    
    # expressions <- sapply(constraints, function(con) {
    #   # Detect operator -- order matters
    #   if (grepl("<=", con)) {
    #     op <- "<="
    #   } else if (grepl(">=", con)) {
    #     op <- ">="
    #   } else if (grepl("<", con)) {
    #     op <- "<"
    #   } else if (grepl(">", con)) {
    #     op <- ">"
    #   } else {
    #     stop(paste("No recognised operator in constraint:", con))
    #   }
    #   
    #   parts <- trimws(strsplit(con, op, fixed = TRUE)[[1]])
    #   lhs_r <- substitute_vars(parts[1])
    #   rhs_r <- substitute_vars(parts[2])
    #   
    #   diff_lr <- paste0("(", lhs_r, ") - (", rhs_r, ")")
    #   diff_rl <- paste0("(", rhs_r, ") - (", lhs_r, ")")
    #   
    #   switch(op,
    #          "<" = paste0(
    #            "if (", diff_lr, " >= 0) {\n",
    #            "  error <- error + (", diff_lr, " + 1e-6)^2 * ", penalty_scale, "\n",
    #            "}"
    #          ),
    #          "<=" = paste0(
    #            "if (", diff_lr, " > 0) {\n",
    #            "  error <- error + (", diff_lr, ")^2 * ", penalty_scale, "\n",
    #            "}"
    #          ),
    #          ">" = paste0(
    #            "if (", diff_rl, " >= 0) {\n",
    #            "  error <- error + (", diff_rl, " + 1e-6)^2 * ", penalty_scale, "\n",
    #            "}"
    #          ),
    #          ">=" = paste0(
    #            "if (", diff_rl, " > 0) {\n",
    #            "  error <- error + (", diff_rl, ")^2 * ", penalty_scale, "\n",
    #            "}"
    #          )
    #   )
    # })
    })
    paste(expressions, collapse = "\n")
  }

  constraints_raw <- conditions
  constraints <- constraints_raw[
    !is.na(constraints_raw) & constraints_raw != ""
  ]
  constraints <- trimws(gsub("\\s+", " ", constraints))
  
  penalty_expression <- build_penalty_expression(
    constraints  = constraints,
    parms_names  = c(names(params_fitted), 
                     names(params_fixed),
                     names(params_functions)),
    states_names = c(names(states_fitted),
                     names(states_fixed),
                     names(states_functions),
                     comp_names)
  )

  
  
  return(list(
    parms          = parms_expression,
    initial_states = initial_states_expression,
    states_params  = states_params_expression,
    sir            = sir_expression,
    ret            = return_expression,
    penalty        = penalty_expression
  ))
}