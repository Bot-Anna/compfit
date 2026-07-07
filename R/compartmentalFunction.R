
compartmentalFunction <- function(modelParams, 
                                  compartment_structure,
                                  return_expression, 
                                  sir_expression,
                                  multipleWorkers = FALSE) {
  # ---- Basic set-up ---- 
  
  ### Initialisation of initial states as well as the coefficients of the first
  ### order terms as well as the coefficients of the second order terms.
  
  # For future reference, here the local names of the compartments X1, ... Xn 
  # are recorded. 
  # IMPORTANT:  The compartments representing susceptible people have to be listed
  #             first!
  #   Inital_values:    initial values as first column
  #   Linearj:          coefficients of linear terms (one column per compartment)
  #   Quadraticj:       coefficients of the quadratic terms containing Xj
  data_vals_coeffs <- modelParams
  
  
  ## ---- Number of compartments ----
  # We determine the number of compartments as they are stored per level
  compartment_cols <- compartment_structure$compartment_cols
  compartment_values <- compartment_structure$compartment_values
  number_of_comps <- compartment_structure$number_of_comps
  comp_names <- compartment_structure$comp_names   # canonical order (States col)

  
  ## ---- Time steps ----
  time_info <- data_vals_coeffs$Others
  time_info <- time_info[time_info != "" & !is.na(time_info)]
  time_info <- gsub(" ", "", time_info)
  
  time_info <- time_info[!is.na(time_info)]
  startpoint <- extract_param(time_info, "startpoint")
  endpoint <- extract_param(time_info, "endpoint")
  partition <- extract_param(time_info, "partition")
  cutoff <- extract_param(time_info, "cutoff")
  # cutoff is OPTIONAL: absent -> Inf ("no cutoff", a harmless literal in the
  # generated R/Julia code that is only referenced by models that use it).
  if (length(cutoff) == 0 || all(is.na(cutoff))) cutoff <- Inf

  time <- seq(0, endpoint-startpoint+1, length.out = ((partition)*(endpoint-startpoint+1)+1))
  time <- as.numeric(time)
  cutoff <- transform_time(cutoff, startpoint)
  # date <- as.Date(seq(as.Date(paste(c(startpoint, "-12-31"), collapse = "")), 
  #                     as.Date(paste(c(endpoint, "-12-31"), collapse = "")), 
  #                     length.out = partition*(endpoint-startpoint)+1))
  
  # date <- as.Date(seq(as.Date(paste(c(startpoint-1, "-12-31"), collapse = "")), 
  #                     as.Date(paste(c(endpoint, "-12-31"), collapse = "")), 
  #                     length.out = partition*(endpoint-startpoint+1)+1))
  
  start_date <- as.Date(paste0(startpoint - 1, "-12-31"))
  end_date   <- as.Date(paste0(endpoint, "-12-31"))
  total_days <- as.numeric(end_date - start_date)
  date       <- start_date + round(seq(0, total_days, length.out = partition*(endpoint-startpoint+1)+1))
  
  ## ---- Variance ----
  # The ones that aren't allowed to have much variance (near constant)
  # penalty_info <- data_vals_coeffs$Variance
  # penalty_info <- penalty_info[penalty_info != "" & !is.na(penalty_info)]
  # 
  # expressions_variance <- sub("^(.*?)\\s*\\(.*\\)$", "\\1", penalty_info)
  # numbers_variance <- as.numeric(sub("^.*\\(((?:\\d*\\.?\\d+)(?:/(?:\\d*\\.?\\d+))?)\\)$", "\\1", penalty_info))

  
  ## ---- Number of levels and their compartments ----
  # A "level" is a subpopulation the mixing/normalisation happens within (age,
  # sex, region, ...), declared by any '_'-prefixed column (the name after '_'
  # is free -- _Level1, _Age1, _Region1 all work; detection is by the '_'
  # prefix, not the word "Level").
  level_compartments <- list()
  if (length(compartment_cols) == 0) {
    # No level column at all: default to a SINGLE level holding every
    # compartment (the standard single-population case), auto-initialised so
    # N1 = total_pop and quadratic terms normalise by the whole population.
    level_compartments <- list(seq_len(number_of_comps))
  } else {
    for (level in seq_along(compartment_cols)) {
      vals <- data_vals_coeffs[[compartment_cols[[level]]]]
      vals <- vals[!is.na(vals) & trimws(as.character(vals)) != ""]
      # Level columns hold compartment references -- a number (1..n) OR a
      # compartment name from the States column. Map both to canonical integer
      # indices via the registry: downstream these are used as numeric membership
      # tests (i %in% vec) and as X[j] subscripts, which must be plain integers.
      level_compartments[[level]] <- .comp_index(vals, comp_names)
    }
  }
  
  ## ---- Model parameters of first order ----
  # First-order coefficients live in ONE column per compartment: Linear1,
  # Linear2, ..., Linear<n> (mirroring the Quadratic<j> layout; replaces the old
  # single column-major `Linear` column). Linear<j> holds the n coefficients of
  # column j of the coefficient matrix; cbind them and transpose so that
  # parameters_first_order[i, j] is the coefficient of X[j] in equation i.
  # A compartment's rate column may be suffixed by its INDEX (Linear1, Linear2,
  # ...) or by its NAME (LinearS, LinearI, ... when States names them S, I, ...).
  # Resolve each compartment i to whichever column actually exists, preferring
  # the numeric form; default to the numeric name.
  .resolve_col <- function(prefix, i) {
    cand <- paste0(prefix, c(i, comp_names[i]))
    hit  <- cand[cand %in% names(data_vals_coeffs)]
    if (length(hit)) hit[1] else cand[1]
  }
  linear_cols <- vapply(seq_len(number_of_comps),
                        function(i) .resolve_col("Linear", i), character(1))
  quad_cols   <- vapply(seq_len(number_of_comps),
                        function(i) .resolve_col("Quadratic", i), character(1))
  # A missing Linear<j>/Quadratic<j> column is a valid way to say "this compartment
  # has no first-/second-order terms": it is treated as all-zero, exactly like a
  # present column full of blanks, and silently (an actual typo in a column name is
  # still caught upstream by validate_modelParams() as an unrecognised column).
  first_order_mat <- vapply(linear_cols, function(cn) {
    col <- if (cn %in% names(data_vals_coeffs)) as.character(data_vals_coeffs[[cn]]) else character(0)
    v <- col[seq_len(number_of_comps)]                                  # first n
    v[is.na(v) | trimws(v) == ""] <- "0"                                # blanks -> 0
    v
  }, character(number_of_comps))                       # n x n, column j = Linear<j>
  parameters_first_order <- t(first_order_mat)
  
  ## ---- Model parameters of second order ----
  # NOTE: The way it should be implemented is very general, however for 
  #       epidemiological models, it is overkill, so we leave it at this easier
  #       version. Due to this annoying detail, the susceptible compartments
  #       have to come first among the Xi's.
  parameters_second_order <- vector("list", number_of_comps)
  for (i in seq_len(number_of_comps)) {
    string <- quad_cols[i]
    # Use only the first n^2 rows (the sheet is rectangular, so a Quadratic column
    # may carry padding rows when another group is longer than n^2). A missing
    # column is treated as all-zero (warned about above).
    col <- if (string %in% names(data_vals_coeffs))
             as.vector(data_vals_coeffs[[string]])[seq_len(number_of_comps^2)]
           else rep("0", number_of_comps^2)
    parameters_second_order[[i]] <- matrix(col, number_of_comps, number_of_comps)
  }
  parameters_second_order <- lapply(parameters_second_order, function(m) {
    # Empty cells mean "no term": null them out. Catch both NA (the usual form
    # from readxl/readr) and "" (possible from a base reader) so a blank is never
    # mistaken for a coefficient. Matches the Linear<j> blank -> 0 handling.
    m[is.na(m) | trimws(as.character(m)) == ""] <- 0
    m
  })
  
  # ---- Equations ----
  ## ---- Build expressions OF LINEAR TERMS one needs for the SIR model ----
  
  # Help vector (initalised as empty) storing all the strings 
  # we build in the upcoming for loop
  vec_help_first_order <- vector()
  
  # For loop that generates a vector containing all the equations needed for the
  # LINEAR part
  for (i in 1:number_of_comps) {
    # Extracts the ith row of the first order term matrix and 
    # joins it with our indicator vector "vec_help_first_order"
    joined_vec <- rbind(parameters_first_order[i,], 1:number_of_comps)
    
    # Remove all terms where the coefficient would be zero.
    joined_vec <- joined_vec[,joined_vec[1,]!=0]
    
    if(length(joined_vec)!=0L){
      # Technical if statement needed, since R automatically converts
      # a column vector into a row vector
      if(NCOL(joined_vec)==1){
        joined_vec <- matrix(joined_vec)
      }
      
      # Helper string
      current_string <- ""
      
      # Loop which appends at each step "coefficient_j*Xj+"
      # As some of the coefficient may be a function, we initialise all the
      # coefficients as a function if the helpstring starts with a leading 
      # asterisk
      for(j in 1:length(joined_vec[1,])){
        help_string <- joined_vec[1,j]
        if(has_leading_dollar_sign(help_string)){
          help_string <- substring(help_string, 2)
          # Build a function containing the coefficient as return expression
          help_function <- function(x){}
          help_function <- funins(help_function, parse(text=help_string),1)
          assign(paste0("f",i,j), help_function)
          # Build string which creates the expression in the equation
          extra_string <- paste(c("f",i,j,"(time)*", 
                                  "X[",joined_vec[2,j],"]+"), collapse="")
        } else {
          extra_string <- paste(c(help_string,"*", 
                                  "X[",joined_vec[2,j],"]+"), collapse="")
        }
        
        current_string <- paste(c(current_string, extra_string), collapse="")
      }
      vec_help_first_order <- append(vec_help_first_order, current_string)
    } else { vec_help_first_order <- append(vec_help_first_order, "") }
  }
  
  ## ---- Build expressions OF QUADRATIC TERMS one needs for the SIR model ----
  
  # Help vector (initalised as empty) storing all the strings 
  # we build in the upcoming for loop
  vec_help_second_order <- vector()
  # To avoid numerical instability, the large second order terms will be initialised
  # at the beginning of the model function, and this term is then added/subtracted 
  # at the appropriate places. To correctly initialise it, we build a vector of
  # these expressions
  vec_help_expressions_second_order <- vector()
  # To understand where the second order terms need to be subtracted from, we 
  # use this additional helping vector add_to_string
  add_to_string <- setNames(numeric(0), character(0))
  # For loop that generates a vector containing all the equations needed for the 
  for (i in 1:number_of_comps) {
    # Determine current level of compartment
    current_level <- which(sapply(level_compartments, 
                                    function(vec) i %in% vec))
    
    # Helper string
    current_string <- ""
    if(as.character(i) %in% names(add_to_string)) {
      current_string <- paste(c(current_string, 
                                add_to_string[as.character(i)]), 
                              collapse="")
    }
    for (j in 1:number_of_comps) {
      # Determines the level of the other compartment
      other_current_level <- which(sapply(level_compartments, 
                                            function(vec) j %in% vec))
      
      # Extracts the ith row of the first order term matrix and 
      # joins it with our indicator vector "vec_help_second_order"
      joined_vec <- rbind(parameters_second_order[[i]][j,], 1:number_of_comps)
      
      # Remove all terms where the coefficient would be zero.
      joined_vec <- joined_vec[,joined_vec[1,]!=0]
      
      if(length(joined_vec)!=0L){
        # Technical if statement needed, since R automatically converts
        # a column vector into a row vector
        if(NCOL(joined_vec)==1){
          joined_vec <- matrix(joined_vec)
        }
        
        # Loop which appends at each step "coefficient_j*Xi*Xj+"
        for(k in 1:length(joined_vec[1,])){
          help_string <- joined_vec[1,k]
          # We need to know to which other compartment to add the expression;
          # this is indicated by *number* at the beginning, which will tell us
          # to which Xi to add it to
          if(has_leading_asterisk(help_string)){
            goes_to <- sub("^\\*(\\d+)\\*.*$", "\\1", help_string)
            help_string <- sub("^\\*\\d+\\*(.*)$", "\\1", help_string)
            
            # The dollar sign indicates if we need to build a function first
            # The functions, for now, are only allowed to have the argument "time"
            if(has_leading_dollar_sign(help_string)){
              help_string <- substring(help_string, 2)
              # Build a function containing the coefficient as return expression
              help_function <- function(x){}
              help_function <- funins(help_function, parse(text=help_string),1)
              assign(paste0("g",i,j), help_function)
              # Build string which creates the expression in the equation
              # Distinguish between both being from same level or both being
              # from different levels 
              if (current_level == other_current_level) {
                extra_string <- paste(c("(g",i,j,"(time)*", 
                                        "X[", j, "]*",
                                        "X[", i, "])/",
                                        "N", current_level), 
                                      collapse="") # instead of i joined_vec[2,k]?
                vec_help_expressions_second_order <- append(
                  vec_help_expressions_second_order,
                  paste(c("secOrd_", i, "_", j, "=", extra_string),
                        collapse=""))
                extra_string <- paste(c("secOrd_", i, "_", j, "+"), collapse="")
              } else {
                extra_string <- paste(c("(g",i,j,"(time)*", 
                                        "X[", j, "]*",
                                        "X[", i,"])/",
                                        "(N", other_current_level, ")"), 
                                      collapse="") # instead of i joined_vec[2,k]?
                vec_help_expressions_second_order <- append(
                  vec_help_expressions_second_order, 
                  paste(c("secOrd_", i, "_", j, "=", extra_string 
                          # "*dampingBoth"), 
                  ),
                  collapse=""))
                extra_string <- paste(c("secOrd_", i, "_", j, "+"), collapse="")
              }
              
            } else {
              if (current_level == other_current_level) {
                extra_string <- paste(c("(", help_string, "*", 
                                        "X[", j, "]*",
                                        "X[", i, "])/",
                                        "N", current_level), 
                                      collapse="")
                vec_help_expressions_second_order <- append(
                  vec_help_expressions_second_order, 
                  paste(c("secOrd_", i, "_", j, "=", extra_string
                          # "*damping", current_level), 
                  ),
                  collapse=""))
                extra_string <- paste(c("secOrd_", i, "_", j, "+"), collapse="")
              } else {
                extra_string <- paste(c("(", help_string, "*", 
                                        "X[", j, "]*",
                                        "X[", i,"])/",
                                        "(N", other_current_level, ")"), 
                                      collapse="")
                vec_help_expressions_second_order <- append(
                  vec_help_expressions_second_order, 
                  paste(c("secOrd_", i, "_", j, "=", extra_string
                          # "*dampingBoth"), 
                  ),
                  collapse=""))
                extra_string <- paste(c("secOrd_", i, "_", j, "+"), collapse="")
              }
              if(is.na(add_to_string[goes_to])){
                current_expr <- ""
              } else {
                current_expr <- add_to_string[goes_to]
              }
              updated_expr <- paste(current_expr, paste0("-", extra_string), sep = "")
              add_to_string[goes_to] <- updated_expr
            }
            goes_to <- ""
            current_string <- paste(c(current_string, extra_string), collapse="")
          }
        }
      } 
    }
    vec_help_second_order <- append(vec_help_second_order, current_string)
  }
  
  ### ---- Outside system ----
  # seq_along (not 1:length): a purely LINEAR model has zero second-order
  # expressions, and 1:length(...) would become 1:0 = c(1, 0), indexing [0] and
  # tripping "argument is of length zero". seq_along() yields integer(0) -> skip.
  for (i in seq_along(vec_help_expressions_second_order)) {
    if (!is.na(vec_help_expressions_second_order[i]) & vec_help_expressions_second_order[i] != "") {
      vec_help_expressions_second_order[i] <- remove_trailing_plus(vec_help_expressions_second_order[i])
    }
  }
  
  
  ## ---- Zeroth-order (constant) terms ----
  # A single optional `Constant` column: row i is a constant term added directly
  # to dX[i] with NO state multiplier -- e.g. a constant inflow / birth /
  # immigration rate, which cannot be expressed as a Linear `*X[j]` or a
  # Quadratic `*X[i]*X[j]` term. Each cell follows the same grammar as a
  # Linear<j> coefficient: a number, a parameter name, an expression of
  # parameters, or a $-prefixed time-function. A blank/missing cell (or an
  # absent column) means 0. Assigned into this function's environment (like the
  # f/g coefficient functions) so the R model closes over it; folded into
  # vec_main below so every backend (R, Julia, Stan) picks it up.
  constant_col <- if ("Constant" %in% names(data_vals_coeffs))
                    as.character(data_vals_coeffs[["Constant"]])[seq_len(number_of_comps)]
                  else rep("0", number_of_comps)
  constant_col[is.na(constant_col) | trimws(constant_col) == ""] <- "0"
  constant_col <- gsub(" ", "", constant_col)
  vec_help_constant <- character(number_of_comps)
  for (i in seq_len(number_of_comps)) {
    val <- constant_col[i]
    if (val == "0") { vec_help_constant[i] <- ""; next }
    if (has_leading_dollar_sign(val)) {
      help_function <- function(x){}
      help_function <- funins(help_function, parse(text = substring(val, 2)), 1)
      assign(paste0("cst", i), help_function)
      vec_help_constant[i] <- paste0("cst", i, "(time)+")
    } else {
      vec_help_constant[i] <- paste0(val, "+")
    }
  }

  ## ---- Builds a combined vector ----

  vec_main <- paste(vec_help_constant, vec_help_first_order, vec_help_second_order)
  for (i in 1:length(vec_main)) {
    vec_main[i] <- remove_trailing_plus(vec_main[i])
  }
  vec_main <- reduce_expression(vec_main)
  # A compartment with NO first- or second-order terms (e.g. every rate column
  # feeding it is blank or absent) has an empty derivative expression. Emit an
  # explicit "0" so the generated code is `dX[i] = 0` rather than a malformed
  # `dX[i] = ` on every backend.
  vec_main[!nzchar(trimws(vec_main))] <- "0"
  
  ## ---- Builds all additional functions for the SIR model ----
  functions_vector <- data_vals_coeffs$Functions
  functions_vector <- functions_vector[functions_vector != "" & 
                                         !is.na(functions_vector)]
  functions_vector <- gsub(" ", "", functions_vector)
  functions_vector <- gsub("<-", "=", functions_vector, fixed = TRUE)
  
  functions_expression <- ""
  
  if (length(functions_vector) != 0) {
    for (i in 1:length(functions_vector)) {
      functions_expression <- paste(c(functions_expression, "\n",
                                      functions_vector[i]), 
                                    collapse = "")
    }
  }
  # For Julia
  functions_expression <- gsub("\\btime\\b", "t", functions_expression)
  
  # ---- SIR ----
  
  ## ---- Builds the SIR model function ----
  ### Using the following expressions:
  ###    "vec_help_first_order"
  ###    "vec_help_second_order"
  ###    "vec_help_expressions_second_order"
  ###    "return_expression"
  ###    "sir_expression"
  
  ## ----Builds master expression ----
  # Containing a call to the function "with", as well as
  # its body filled with the expressions, and "return_expression"
  master_expression <- ""
  for (i in seq_along(level_compartments)) {
    master_expression <- paste(c(master_expression,"\n N", i, "= ("), collapse="")
    for (j in level_compartments[[i]]) {

      master_expression <- paste(c(master_expression,"X[", j, "]+"), collapse="")
    }
    master_expression <- remove_trailing_plus(master_expression)
    master_expression <- paste(c(master_expression, ")"),
                               collapse = "")
  }

  master_expression <- paste(c(master_expression,
                               "\n total_pop = "),
                             collapse = "")
  for (i in seq_along(level_compartments)) {
    master_expression <- paste(c(master_expression,
                                 "N", i, "+"),
                               collapse = "")
  }
  master_expression <- remove_trailing_plus(master_expression)
  
  # Adds the expression so it's compatible with Julia
  master_expression <- paste(c(master_expression, "\n", sir_expression),
                             collapse = "")
  
  # Adds all the functions that are set
  master_expression <- paste(c(master_expression, functions_expression), 
                             collapse = "")
  
  # Adds the second order terms using vec_help_expressions_second_order.
  # seq_along (not 1:length): empty for a purely linear model -> skip cleanly.
  for (i in seq_along(vec_help_expressions_second_order)) {
    master_expression <- paste(c(master_expression,
                                 "\n ",
                                 vec_help_expressions_second_order[[i]]),
                               collapse="")
  }
  
  # Adds all the differential equations
  for (i in 1:number_of_comps) {
    master_expression <- paste(c(master_expression, paste(c("\ndX",i,"=",
                                                            vec_main[i]),
                                                          collapse = "")),
                               collapse = "")
  }
  
  master_expression <- remove_trailing_plus(master_expression)
  
  # Adds the return expression at the end
  master_expression <- paste(c(master_expression, return_expression), 
                             collapse = "")
  
  
  
  # Empty SIR model function which is filled in the next step using funins
  sir_model <- function(X, p, t){
  }
  sir_model <- funins(sir_model, parse(text=master_expression), 1)
  
  # Build the Julia version of the ODE function
  julia_code <- buildJuliaODEFunction(
    sir_expression = sir_expression,
    functions_expression = functions_expression,
    vec_help_expressions_second_order = vec_help_expressions_second_order,
    vec_main = vec_main,
    number_of_comps = number_of_comps,
    level_compartments = level_compartments,
    cutoff = cutoff,
    startpoint = startpoint
    )

  # Build the Stan version of the ODE (same IR). Wrapped so a codegen hiccup on
  # an exotic model never breaks the MLE / Julia / R paths -- stan_code is NULL
  # then, and the Stan backend reports it clearly when actually requested.
  stan_code <- tryCatch(
    buildStanODEFunction(
      sir_expression = sir_expression,
      functions_expression = functions_expression,
      vec_help_expressions_second_order = vec_help_expressions_second_order,
      vec_main = vec_main,
      number_of_comps = number_of_comps,
      level_compartments = level_compartments,
      cutoff = cutoff,
      startpoint = startpoint
    ), error = function(e) NULL)

  # ---- Return statement ----
  # We return both the function and the vector of functions needed for multiple
  # workers, or just the former, depending on the value of multipleWorkers
  if(multipleWorkers == TRUE) {
      return(list(compartmental_function = sir_model,
                  date                   = date,
                  julia_code             = julia_code,
                  stan_code              = stan_code,
                  vector_of_functions    = vector_of_functions))
    } else {
      return(list(compartmental_function = sir_model,
                  date                   = date,
                  julia_code             = julia_code,
                  stan_code              = stan_code))
    }
  
}




## ============================================================
## Julia ODE function builder
## Drop-in replacement for the R callback passed to de$ODEProblem
## 
## USAGE:
##   julia_func <- buildJuliaODEFunction(
##                   sir_expression,         # from generateExpressions()
##                   functions_expression,   # built inside compartmentalFunction()
##                   vec_help_expressions_second_order,
##                   vec_main,
##                   number_of_comps,
##                   level_compartments,
##                   cutoff,
##                   startpoint
##                 )
##   julia_eval(julia_func)   # defines compartmental_function_jl in Julia
##
##   prob <- julia_eval("ODEProblem(compartmental_function_jl, X, t, p)")
##   sol  <- de$solve(prob, de$BS3(), saveat = time_integer,
##                    abstol = 1e-4, reltol = 1e-4)
## ============================================================


buildJuliaODEFunction <- function(sir_expression,
                                   functions_expression,
                                   vec_help_expressions_second_order,
                                   vec_main,
                                   number_of_comps,
                                   level_compartments,
                                   cutoff,
                                   startpoint) {

  ## ----------------------------------------------------------
  ## 1. Helper: convert R expression strings to Julia syntax
  ## ----------------------------------------------------------
  r_to_julia <- function(s) {
    # R uses ^ for power, Julia uses ^  (same — no change needed)
    # R uses if_else(cond, a, b) — Julia uses cond ? a : b
    # We do a simple recursive substitution for the nested if_else
    # pattern that appears in your time-varying lookup tables.
    
    # if_else(cond, a, b)  ->  (cond ? a : b)
    # Applied repeatedly to handle nesting
    max_iter <- 200
    iter <- 0
    while (grepl("if_else\\(", s) && iter < max_iter) {
      s <- gsub_if_else(s)
      iter <- iter + 1
    }
    
    # ifelse(cond, a, b) — same treatment
    iter <- 0
    while (grepl("\\bifelse\\(", s) && iter < max_iter) {
      s <- gsub_ifelse(s)
      iter <- iter + 1
    }
    
    # TRUE/FALSE -> true/false
    s <- gsub("\\bTRUE\\b",  "true",  s)
    s <- gsub("\\bFALSE\\b", "false", s)
    
    # R: <- is assignment, Julia uses =
    s <- gsub("<-", "=", s, fixed = TRUE)
    
    s
  }
  
  ## Replaces the OUTERMOST if_else(...) with Julia ternary
  gsub_if_else <- function(s) {
    # Find "if_else(" and walk to matching close paren, 
    # then split into 3 args
    pattern <- "if_else\\("
    m <- regexpr(pattern, s)
    if (m == -1) return(s)
    
    start <- m[1] + attr(m, "match.length")[1] - 1  # position of opening (
    args <- extract_three_args(s, start)
    if (is.null(args)) return(s)  # bail if parse fails
    
    replacement <- paste0("(", args[1], " ? ", args[2], " : ", args[3], ")")
    
    before <- substr(s, 1, m[1] - 1)
    after  <- substr(s, args[[4]], nchar(s))  # after closing )
    paste0(before, replacement, after)
  }
  
  gsub_ifelse <- function(s) {
    pattern <- "\\bifelse\\("
    m <- regexpr(pattern, s)
    if (m == -1) return(s)
    
    start <- m[1] + attr(m, "match.length")[1] - 1
    args <- extract_three_args(s, start)
    if (is.null(args)) return(s)
    
    replacement <- paste0("(", args[1], " ? ", args[2], " : ", args[3], ")")
    before <- substr(s, 1, m[1] - 1)
    after  <- substr(s, args[[4]], nchar(s))
    paste0(before, replacement, after)
  }
  
  ## Extracts 3 comma-separated arguments from a parenthesised expression.
  ## 'pos' is the position of the opening parenthesis in 's'.
  ## Returns list(arg1, arg2, arg3, pos_after_close).
  extract_three_args <- function(s, pos) {
    chars  <- strsplit(s, "")[[1]]
    depth  <- 0
    args   <- character(3)
    arg_i  <- 1
    buf    <- ""
    i      <- pos  # points at '('
    
    if (chars[i] != "(") return(NULL)
    i <- i + 1  # step past '('
    depth <- 1
    
    while (i <= length(chars) && arg_i <= 3) {
      ch <- chars[i]
      if (ch == "(") {
        depth <- depth + 1
        buf <- paste0(buf, ch)
      } else if (ch == ")") {
        depth <- depth - 1
        if (depth == 0) {
          args[arg_i] <- trimws(buf)
          return(list(args[1], args[2], args[3], i + 1))
        }
        buf <- paste0(buf, ch)
      } else if (ch == "," && depth == 1) {
        args[arg_i] <- trimws(buf)
        arg_i <- arg_i + 1
        buf <- ""
      } else {
        buf <- paste0(buf, ch)
      }
      i <- i + 1
    }
    NULL  # parse failed
  }

  ## ----------------------------------------------------------
  ## 2. Build N-population sum lines  (N1 = X[1]+...+X[k])
  ## ----------------------------------------------------------
  n_lines <- ""
  for (i in seq_along(level_compartments)) {
    idx <- level_compartments[[i]]
    terms <- paste0("X[", idx, "]", collapse = "+")
    n_lines <- paste0(n_lines, "\n    N", i, " = ", terms)
  }
  all_N <- paste0("N", seq_along(level_compartments), collapse = "+")
  n_lines <- paste0(n_lines, "\n    total_pop = ", all_N)

  ## ----------------------------------------------------------
  ## 3. sir_expression — parameter unpacking (already built by
  ##    generateExpressions as "param = p[i]" lines)
  ## ----------------------------------------------------------
  sir_expr_julia <- r_to_julia(sir_expression)

  ## ----------------------------------------------------------
  ## 4. functions_expression — time-varying lookup tables
  ##    (h_nCH, perc_immH, etc.)
  ## ----------------------------------------------------------
  func_expr_julia <- r_to_julia(functions_expression)
  
  ## Inject cutoff and startpoint as literals (they are R-side 
  ## constants that the Julia function needs to see)
  func_expr_julia <- gsub("\\bcutoff\\b",    as.character(cutoff),    func_expr_julia)
  func_expr_julia <- gsub("\\bstartpoint\\b", as.character(startpoint), func_expr_julia)

  ## ----------------------------------------------------------
  ## 5. Second-order (secOrd) initialisation expressions
  ## ----------------------------------------------------------
  sec_ord_lines <- ""
  for (expr in vec_help_expressions_second_order) {
    if (!is.na(expr) && nchar(trimws(expr)) > 0) {
      sec_ord_lines <- paste0(sec_ord_lines, "\n    ", r_to_julia(expr))
    }
  }

  ## ----------------------------------------------------------
  ## 6. Differential equations  dX1 = ...
  ## ----------------------------------------------------------
  dX_lines <- ""
  for (i in seq_len(number_of_comps)) {
    dX_lines <- paste0(dX_lines, "\n    dX[", i, "] = ", r_to_julia(vec_main[i]))
  }

  ## ----------------------------------------------------------
  ## 7. Assemble complete Julia function string
  ## ----------------------------------------------------------
  julia_code <- paste0(
    'function compartmental_function_jl(dX, X, p, t)\n',
    n_lines, "\n",
    sir_expr_julia, "\n",
    func_expr_julia, "\n",
    sec_ord_lines, "\n",
    dX_lines, "\n",
    '    return nothing\n',
    'end'
  )

  return(julia_code)
}


## ============================================================
## Wrapper: replaces the de$ODEProblem(compartmental_function, ...)
## call in lossFunction.R with a pure-Julia version.
##
## Call this ONCE after compartmentalFunction() returns, before
## constructing the loss function.
## ============================================================

#' Register the ODE function in Julia
#'
#' Low-level solver bridge: defines the generated pure-Julia ODE function
#' (`compartmental_function_jl`) in the running Julia session so
#' [solveWithJulia()] can call it. Exported because the self-contained scripts
#' emitted by [extract_code()] call it after `library(compfit)`.
#'
#' @param julia_code Character; the Julia ODE source (e.g. `get_julia_code(fit)`).
#' @return `TRUE`, invisibly (called for its side effect in the Julia session).
#' @export
registerJuliaODEFunction <- function(julia_code) {
  .compfit_ensure_julia()
  # Make OrdinaryDiffEq available (already loaded by diffeqr but
  # explicit import makes solve/ODEProblem available at top level)
  JuliaCall::julia_command("using OrdinaryDiffEq")
  JuliaCall::julia_eval(julia_code)
  message("Julia ODE function 'compartmental_function_jl' registered.")
}


## ============================================================
## Replacement solve call for lossFunction.R
##
## Replace:
##   prob <- de$ODEProblem(compartmental_function, X, t, p)
##   sol  <- de$solve(prob, de$BS3(), saveat = time_integer, ...)
##
## With:
##   sol <- solveWithJulia(X, t, p, time_integer)
## ============================================================

#' Solve the ODE in Julia
#'
#' Low-level solver bridge: solves the registered Julia ODE at a given initial
#' state, parameter vector and save times. Mirrors [solveWithR()] (the deSolve
#' backend) and returns the same shape. Exported because the scripts emitted by
#' [extract_code()] call it.
#'
#' @param X Named numeric initial-state vector.
#' @param t Length-2 numeric integration span `c(t0, t1)`.
#' @param p Numeric parameter vector (in ODE order).
#' @param time_integer Numeric vector of save times.
#' @param solver Julia solver expression.
#' @param abstol,reltol Solver tolerances.
#' @return A list with `matrix` (states by time) and `t` (the save times).
#' @export
solveWithJulia <- function(X, t, p, time_integer,
                            solver  = "BS3()",
                            abstol  = 1e-8,
                            reltol  = 1e-8) {
  .compfit_ensure_julia()

  JuliaCall::julia_assign("_X_jl",   X)
  JuliaCall::julia_assign("_t_jl",   t)
  JuliaCall::julia_assign("_p_jl",   p)
  JuliaCall::julia_assign("_save_jl", time_integer)
  
  cmd <- paste0(
    "_prob_jl = ODEProblem(compartmental_function_jl, _X_jl, (_t_jl[1], _t_jl[2]), _p_jl); ",
    "_sol_jl = solve(_prob_jl, ", solver, ", ",
    "saveat=_save_jl, abstol=", abstol, ", reltol=", reltol, "); nothing"
  )
  invisible(capture.output(
    JuliaCall::julia_command(cmd),
    type = "output"
  ))
  
  
  # cmd <- paste0(
  #   "_prob_jl = ODEProblem(compartmental_function_jl, _X_jl, (_t_jl[1], _t_jl[2]), _p_jl); ",
  #   "_sol_jl  = solve(_prob_jl, ", solver, ", ",
  #   "saveat=_save_jl, abstol=", abstol, ", reltol=", reltol, ")"
  # )
  # # JuliaCall::julia_command(cmd)
  # 
  # invisible(JuliaCall::julia_command(cmd))
  
  sol_matrix <- JuliaCall::julia_eval("Matrix(_sol_jl)")
  sol_t      <- JuliaCall::julia_eval("_sol_jl.t")

  list(matrix = sol_matrix, t = sol_t)
}

## ===================== PURE-R SOLVER BACKEND =====================
## A Julia-free alternative to registerJuliaODEFunction()/solveWithJulia(),
## using the R compartmental_function closure + deSolve. Mirrors the Julia
## pair's interface so the loss closure and solve_and_evaluate() can call
## either backend interchangeably. Slower than Julia, but needs no Julia
## runtime -- the basis of solver_control(backend = "r").

## Register the R ODE closure for the current session (mirrors
## registerJuliaODEFunction). solveWithR() reads it from the package state env.
registerRODEFunction <- function(compartmental_function) {
  .compfit_state$r_ode <- compartmental_function
  invisible(TRUE)
}

## Drop-in replacement for solveWithJulia() using deSolve::ode (lsoda: an
## auto-switching stiff/non-stiff method, a good analogue of the Julia default).
## Returns list(matrix = states x time, t = save points), exactly like
## solveWithJulia(), so downstream code is agnostic to the backend.
# deSolve integration methods; used to validate a solver_control(solver=) string
# on the R backend. Anything else (e.g. the default Julia solver expression)
# falls back to the auto-switching `lsoda`.
.DESOLVE_METHODS <- c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
  "euler", "rk4", "ode23", "ode45", "radau", "bdf", "bdf_d", "adams",
  "impAdams", "impAdams_d", "iteration")
.desolve_method <- function(solver_str) {
  s <- as.character(solver_str)[1]
  if (length(s) && !is.na(s) && s %in% .DESOLVE_METHODS) s else "lsoda"
}

solveWithR <- function(X, t, p, time_integer, abstol = 1e-8, reltol = 1e-8,
                       method = "lsoda") {
  cf <- .compfit_state$r_ode
  if (is.null(cf))
    stop("No R ODE function registered. Build/fit with solver_control(backend = \"r\"), ",
         "or load a fit saved from an R-backend run.")
  if (!requireNamespace("deSolve", quietly = TRUE))
    stop("solveWithR() needs the 'deSolve' package (install.packages(\"deSolve\")).")

  deriv <- function(t, y, parms) list(cf(y, parms, t))   # deSolve wants list(dY)
  out <- deSolve::ode(y = X, times = as.numeric(time_integer), func = deriv,
                      parms = p, method = method, atol = abstol, rtol = reltol)

  states <- t(out[, -1, drop = FALSE])   # times x (1+nstate) -> states x time
  rownames(states) <- names(X)
  list(matrix = states, t = out[, 1])
}
