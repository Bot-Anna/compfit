
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
  # sex, region, ...), declared by a `Level_<id>` column (the id after `Level_`
  # is free -- Level_1, Level_LA, Level_Age1 all work; it names the level and
  # drives the companions N_<id> / Mixing_<id> / Pool_<id>).
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

  ## ---- Per-level mixing pools (full-pool transmission denominators) ----
  # A level column `_<name>` may carry a companion that REPLACES the raw headcount
  # N<L> in the transmission (second-order) denominators for that level with a
  # weighted / custom "mixing pool" Nw<L>. Two mutually-exclusive forms:
  #   * `Mixing_<name>` -- one weight per compartment, aligned to the `Level_<name>`
  #     column BY ROW: the weight in row i is for the compartment named in row i of
  #     Level_<name>. A blank weight on an in-level compartment defaults to 1, so an
  #     absent / all-blank column reproduces N<L> exactly. A weight in a row with NO
  #     Level_<name> entry pulls the States-position compartment X[row] into the pool
  #     (commuting / contact-matrix mixing). Robust to the States ordering (grouped
  #     OR interleaved); the older row=X[k] alignment mis-weighted an interleaved order.
  #   * `Pool_<name>` -- ONE expression giving the whole pool as a function of the
  #     level head counts (N1..Nk, total_pop), parameters, Functions, and time,
  #     e.g. `N1 + c*N2/(1+N2/K)` (a saturating cross-level pool). Floored at a
  #     small positive value so a pool that momentarily hits 0 cannot divide by
  #     zero; a genuinely non-positive pool is a modelling error (see validate).
  # Raw N<L> and total_pop are left untouched (still the true head counts, and
  # still user-referenceable). Weights/expressions follow the coefficient grammar.
  .NW_FLOOR        <- "1e-8"
  level_names      <- if (length(compartment_cols)) sub("^Level_", "", compartment_cols, ignore.case = TRUE) else character(0)
  # Named aliases for the level head counts: alongside the positional N1..Nk, emit
  # `N_<levelname> = N<i>` (e.g. `N_LA = N1`) so a coefficient / Function / Pool can
  # refer to a stratum by name instead of tracking its column order. Only when the
  # levels are explicitly named (a `_<name>` column exists).
  n_alias_defs <- if (length(compartment_cols))
                    paste0("N_", level_names, " = N", seq_along(level_names))
                  else character(0)
  mixing_cols      <- if (length(level_names)) paste0("Mixing_", level_names) else character(0)
  pool_cols        <- if (length(level_names)) paste0("Pool_",   level_names) else character(0)
  has_mixing       <- mixing_cols %in% names(data_vals_coeffs)
  has_pool         <- pool_cols   %in% names(data_vals_coeffs)
  level_has_mixing <- has_mixing | has_pool              # "has a custom denominator"
  # The denominator token for level L: the custom pool if it has one, else the raw
  # headcount (so models without a Mixing/Pool column are byte-identical).
  denom <- function(L) {
    if (length(level_has_mixing) >= L && level_has_mixing[L]) paste0("Nw", L)
    else paste0("N", L)
  }
  # R-syntax `Nw<L> = ...` definitions, emitted AFTER parameter unpacking and the
  # Functions block by every backend (a weight/expression may use a fitted
  # parameter); N1..Nk and total_pop are already in scope from the top of the body.
  nw_defs <- character(0)
  for (L in seq_along(level_compartments)) {
    if (length(level_has_mixing) < L || !level_has_mixing[L]) next
    if (has_pool[L]) {                                   # Pool_<name>: a pool expression
      pv   <- as.character(data_vals_coeffs[[pool_cols[L]]])
      pv   <- pv[!is.na(pv) & nzchar(trimws(pv))]
      expr <- .names_to_slots(gsub(" ", "", trimws(pv[1])), comp_names)
      rhs  <- paste0("max((", expr, "), ", .NW_FLOOR, ")")   # positive floor guards /0
    } else {                                             # Mixing_<name>: weights aligned to Level_<name>
      # Each weight sits in the SAME ROW as its compartment's entry in the
      # Level_<name> column (how the column is naturally written), so a blank
      # in-level weight defaults to 1 and the pool is robust to the States
      # ordering (grouped OR interleaved). A weight in a row that has NO Level
      # entry pulls the States-position compartment X[row] into the pool
      # (cross-level contact mixing). The old row=X[k] alignment silently
      # mis-weighted an interleaved (non-level-grouped) States order.
      lvlcol <- as.character(data_vals_coeffs[[compartment_cols[[L]]]])
      mixcol <- as.character(data_vals_coeffs[[mixing_cols[L]]])
      terms  <- character(0)
      for (i in seq_len(max(length(lvlcol), length(mixcol)))) {
        lv <- if (i <= length(lvlcol)) trimws(lvlcol[i]) else NA_character_
        mw <- if (i <= length(mixcol)) trimws(mixcol[i]) else NA_character_
        if (!is.na(lv) && nzchar(lv)) {                 # row lists an in-level compartment
          k <- .comp_index(lv, comp_names)
          w <- if (is.na(mw) || !nzchar(mw)) "1" else mw  # blank in-level weight -> 1
        } else if (!is.na(mw) && nzchar(mw)) {          # weight, no Level entry -> pull in X[row]
          k <- i; w <- mw
        } else next
        if (is.na(k) || k < 1L || k > number_of_comps) next
        w <- .names_to_slots(gsub(" ", "", w), comp_names)
        if (w == "0") next
        terms <- c(terms, if (w == "1") paste0("X[", k, "]")
                          else paste0("(", w, ")*X[", k, "]"))
      }
      # All-zero pool would divide by zero; fall back to the raw headcount.
      rhs <- if (length(terms)) paste(terms, collapse = "+") else paste0("N", L)
    }
    nw_defs <- c(nw_defs, paste0("Nw", L, " = ", rhs))
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
      
      # Loop which appends at each step "coefficient_j*Xj+". A coefficient is a
      # number, a parameter name, or an expression of parameters/time/STATES; a
      # time-varying coefficient is written as a named `Functions` entry and
      # referenced here by name. State names in the coefficient are rewritten to
      # their X[k] slot.
      for(j in 1:length(joined_vec[1,])){
        help_string <- .names_to_slots(joined_vec[1,j], comp_names)
        extra_string <- paste(c(help_string,"*",
                                "X[",joined_vec[2,j],"]+"), collapse="")
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
          # this is indicated by *target* at the beginning. The target is a
          # compartment INDEX (*5*) or NAME (*I*, when States name them S/I/R);
          # .comp_index() resolves either to the canonical integer index that
          # keys add_to_string below.
          if(has_leading_asterisk(help_string)){
            goto_token <- sub("^\\*([^*]+)\\*.*$", "\\1", help_string)
            goes_to <- as.character(.comp_index(goto_token, comp_names))
            help_string <- sub("^\\*[^*]+\\*(.*)$", "\\1", help_string)
            # State names in the coefficient -> X[k] (the goto target above is
            # resolved by name separately, so it is left untouched).
            help_string <- .names_to_slots(help_string, comp_names)

            # Same vs different mixing level decides the normalising denominator.
            # denom() yields Nw<level> when that level has a Mixing_<name> column
            # (full-pool weighting) and the raw N<level> headcount otherwise. A
            # time-varying coefficient is a named `Functions` entry referenced here.
            if (current_level == other_current_level) {
              extra_string <- paste(c("(", help_string, "*",
                                      "X[", j, "]*",
                                      "X[", i, "])/",
                                      denom(current_level)),
                                    collapse="")
            } else {
              extra_string <- paste(c("(", help_string, "*",
                                      "X[", j, "]*",
                                      "X[", i,"])/",
                                      "(", denom(other_current_level), ")"),
                                    collapse="")
            }
            vec_help_expressions_second_order <- append(
              vec_help_expressions_second_order,
              paste(c("secOrd_", i, "_", j, "=", extra_string), collapse=""))
            extra_string <- paste(c("secOrd_", i, "_", j, "+"), collapse="")
            if(is.na(add_to_string[goes_to])){
              current_expr <- ""
            } else {
              current_expr <- add_to_string[goes_to]
            }
            updated_expr <- paste(current_expr, paste0("-", extra_string), sep = "")
            add_to_string[goes_to] <- updated_expr
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
  # Linear<j> coefficient: a number, a parameter name, or an expression of
  # parameters (a time-varying constant is a named `Functions` entry referenced
  # here). A blank/missing cell (or an absent column) means 0. Folded into
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
    vec_help_constant[i] <- paste0(.names_to_slots(val, comp_names), "+")
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
  # Sorted into dependency order first: the block below is emitted as straight-
  # line code, so an entry must follow everything it references. Row order in the
  # sheet carries no meaning -- see .function_order(). The SAME ordered vector
  # feeds `derived_spec$functions_raw` further down, so the generated ODE and the
  # replay in .derived_columns() agree.
  functions_all <- as.character(data_vals_coeffs$Functions)
  functions_all <- functions_all[!is.na(functions_all) & nzchar(trimws(functions_all))]
  functions_all <- functions_all[.order_functions(functions_all)]

  functions_vector <- functions_all
  functions_vector <- gsub(" ", "", functions_vector)
  functions_vector <- gsub("<-", "=", functions_vector, fixed = TRUE)
  # Rewrite STATE NAMES to X[k] in each entry's RHS (the LHS is the new variable's
  # own name and is left as-is), so a Functions definition can be written in terms
  # of the compartments, e.g. `Hinf = tau*R_HA + C_HA + delta*D_HA`.
  functions_vector <- vapply(functions_vector, function(e) {
    pos <- regexpr("=", e, fixed = TRUE)
    if (pos < 1L) return(e)
    paste0(substr(e, 1L, pos), .names_to_slots(substring(e, pos + 1L), comp_names))
  }, character(1), USE.NAMES = FALSE)

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

  # Named level head-count aliases (N_<level> = N<i>), right after the N block.
  for (na in n_alias_defs)
    master_expression <- paste(c(master_expression, "\n ", na), collapse = "")

  # Adds the expression so it's compatible with Julia
  master_expression <- paste(c(master_expression, "\n", sir_expression),
                             collapse = "")
  
  # Adds all the functions that are set
  master_expression <- paste(c(master_expression, functions_expression),
                             collapse = "")

  # Adds the weighted-mixing pools Nw<L> (full-pool transmission denominators).
  # After parameter unpacking and the functions block so a weight may itself be a
  # fitted parameter or a Functions entry; before the second-order terms, which
  # divide by them. Empty for a model with no Mixing_<name> column.
  for (nd in nw_defs)
    master_expression <- paste(c(master_expression, "\n ", nd), collapse = "")

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
    nw_defs = nw_defs,
    n_alias_defs = n_alias_defs,
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
      nw_defs = nw_defs,
      n_alias_defs = n_alias_defs,
      cutoff = cutoff,
      startpoint = startpoint
    ), error = function(e) NULL)

  # ---- Derived-quantity replay spec ----
  # Enough to recompute every Functions-column quantity (q_CH, omega_*, p_*,
  # sigmoids, ...) and the level head counts on an ALREADY-SOLVED trajectory, so
  # data/dummy formulas can reference them by name (see solve_and_evaluate()).
  # We keep the RAW Functions (original state names + `time`, not the X[k]/`t`
  # codegen form) because they are replayed as plain R, vectorised over the grid,
  # against the named state columns. `level_members`/`level_names`/`startpoint`/
  # `cutoff` mirror exactly what the generated ODE uses, so the replay agrees with
  # the solve.
  functions_raw <- functions_all          # already in dependency order (see above)
  derived_spec <- list(
    functions_raw = functions_raw,
    level_members = level_compartments,   # list of integer index vectors
    level_names   = level_names,          # character (empty when no Level_ column)
    comp_names    = comp_names,           # canonical compartment order
    startpoint    = startpoint,
    cutoff        = cutoff
  )

  # ---- Return statement ----
  # We return both the function and the vector of functions needed for multiple
  # workers, or just the former, depending on the value of multipleWorkers
  if(multipleWorkers == TRUE) {
      return(list(compartmental_function = sir_model,
                  date                   = date,
                  julia_code             = julia_code,
                  stan_code              = stan_code,
                  derived_spec           = derived_spec,
                  vector_of_functions    = vector_of_functions))
    } else {
      return(list(compartmental_function = sir_model,
                  date                   = date,
                  julia_code             = julia_code,
                  stan_code              = stan_code,
                  derived_spec           = derived_spec))
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
                                   nw_defs = character(0),
                                   n_alias_defs = character(0),
                                   cutoff,
                                   startpoint) {

  ## ----------------------------------------------------------
  ## 1. Helper: convert R expression strings to Julia syntax
  ## ----------------------------------------------------------

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
  # Named level head-count aliases (N_<level> = N<i>).
  for (na in n_alias_defs) n_lines <- paste0(n_lines, "\n    ", na)

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
  ## 4b. Weighted-mixing pools Nw<L> (full-pool denominators)
  ##     After param unpacking + functions, before secOrd uses them.
  ## ----------------------------------------------------------
  nw_lines <- ""
  for (nd in nw_defs) {
    if (!is.na(nd) && nzchar(trimws(nd)))
      nw_lines <- paste0(nw_lines, "\n    ", r_to_julia(nd))
  }

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
    nw_lines, "\n",
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
