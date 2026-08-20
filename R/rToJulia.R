# ============================================================
# rToJulia.R
# R-expression -> Julia-syntax conversion helpers.
#
# These were nested inside buildJuliaODEFunction(); they are package-level here
# so the Bayesian (Turing) code generator can reuse the SAME translation the ODE
# codegen uses, rather than duplicating it. Bodies are unchanged, and the names
# are unchanged, so buildJuliaODEFunction()'s call sites resolve to these via the
# package namespace exactly as before.
# ============================================================

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
