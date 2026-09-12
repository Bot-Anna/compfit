test_that("test-unit-function-order", {
# ============================================================
# test-unit-function-order.R   (pure R; no Julia)
# The Functions column is emitted as a straight-line assignment block, so each
# entry must follow the ones it references. `.function_order()` recovers that
# order from the reference graph instead of trusting sheet row order, and
# reports the three faults no reordering can fix: a duplicate definition, a
# self-reference, and a cycle.
# ============================================================
th_load_pure(c("utils.R", "validate.R"))

nm_of <- function(fr, ord) sub("\\s*<-.*$", "", fr[ord])

th_section("a forward reference is reordered away")
# `Pi_A` written last but used by the first entry -- the exact shape that raised
# UndefVarError in Julia, and that silently dropped the entry in the replay.
fr  <- c("phiLow_A <- (1-Pi_A)/((1-Pi_A)+Pi_A)",
         "sigma_A  <- Pi_A*total_switch",
         "Pi_A     <- Pi_CH")
iss <- .function_order(fr)
chk("no fault is reported", !length(c(iss$duplicates, iss$self, iss$cycle)))
chk_equal("the definition is hoisted above both of its users",
          nm_of(fr, iss$order), c("Pi_A", "phiLow_A", "sigma_A"))

th_section("an already-valid order is returned unchanged")
# Stability matters beyond tidiness: every sheet that works today must keep
# generating byte-identical code, including the `$code` recipe.
fr <- c("b <- 1", "a <- b", "c <- 1", "d <- a + c")
chk_equal("the permutation is the identity", .function_order(fr)$order, 1:4)

th_section("chains and diamonds resolve")
fr  <- c("d <- b + c", "b <- a", "c <- a", "a <- 1")
ord <- .function_order(fr)$order
pos <- setNames(seq_along(ord), nm_of(fr, ord))
chk("every entry follows the ones it references",
    pos[["a"]] < pos[["b"]] && pos[["a"]] < pos[["c"]] &&
      pos[["b"]] < pos[["d"]] && pos[["c"]] < pos[["d"]])

th_section("external symbols constrain nothing")
# Parameters, states, N<level> and the grid symbols are assigned outside the
# block, so referencing them must not create an edge (or a phantom cycle).
fr <- c("omega <- n_H*(N_HighCH + chi*S_P)", "q <- ifelse(time+startpoint < 2015, a, b)")
chk_equal("order is untouched", .function_order(fr)$order, 1:2)

th_section("a duplicate definition is rejected")
iss <- .function_order(c("a <- 1", "b <- a", "a <- 2"))
chk_equal("the duplicated name is named", iss$duplicates, "a")
chk("the message says it is defined more than once",
    grepl("defined more than once", .function_order_errors(iss)[1]))

th_section("a self-reference is rejected")
iss <- .function_order(c("x <- x + 1"))
chk_equal("the self-referencing name is named", iss$self, "x")
chk("the message says a Function is a definition, not an update",
    grepl("not an update", .function_order_errors(iss)[1]))
chk("a self-reference is still caught inside a call",
    identical(.function_order(c("y <- 2", "x <- exp(x)/y"))$self, "x"))

th_section("a cycle is rejected, and the message names the loop")
iss <- .function_order(c("a <- b", "b <- c", "c <- a"))
chk("a cycle is reported", length(iss$cycle) == 3L)
chk("all three members are named", setequal(iss$cycle, c("a", "b", "c")))
msg <- .function_order_errors(iss)[1]
chk("the message says the definition is circular", grepl("circular definition", msg))
chk("it spells out the loop closing back on its start",
    grepl(paste(sQuote(c(iss$cycle, iss$cycle[1])), collapse = " -> "), msg, fixed = TRUE))

th_section("entries merely downstream of a cycle are not blamed for it")
# `d` cannot be placed either, but it is a victim, not a member -- naming it
# would send the user to the wrong row.
iss <- .function_order(c("a <- b", "b <- a", "d <- a"))
chk("only the two cycle members are named", setequal(iss$cycle, c("a", "b")))

th_section("validate_modelParams surfaces the faults")
mk <- function(functions) {
  cols <- list(`Level_1` = c("1", "2"),
               Others = c("startpoint=2000", "endpoint=2005", "partition=4"),
               States = c("*X1=990", "*X2=10"), Functions = functions,
               Parameters = c("beta=[0,2]", "gamma=[0,1]"), Conditions = "",
               Linear1 = c("0", "0"), Quadratic1 = c("0", "*2*-beta"),
               Linear2 = c("0", "-gamma"), Quadratic2 = c("0", "0"))
  nr   <- max(lengths(cols))
  cols <- lapply(cols, function(x) c(as.character(x), rep("", nr - length(x))))
  do.call(data.frame, c(cols, check.names = FALSE, stringsAsFactors = FALSE))
}
errmsg <- function(mp) tryCatch({ validate_modelParams(mp); NA_character_ },
                                error = function(e) conditionMessage(e))

chk_ok("a sheet written in dependency order validates",
       validate_modelParams(mk(c("f<-beta", "g<-f*2"))))
chk_ok("so does the same sheet written out of order -- reordering fixes it",
       validate_modelParams(mk(c("g<-f*2", "f<-beta"))))
chk("a cycle fails validation",
    grepl("circular definition", errmsg(mk(c("f<-g", "g<-f")))))
chk("a duplicate definition fails validation",
    grepl("defined more than once", errmsg(mk(c("f<-beta", "f<-gamma")))))
chk("a self-reference fails validation",
    grepl("not an update", errmsg(mk(c("f<-f+1", "")))))

th_section("row order cannot change what the model computes")
# The real guarantee: permuting the Functions column leaves the generated ODE
# numerically identical. (A topological order is not unique, so the generated
# TEXT may differ between permutations -- only the computed values must not.)
fn <- c("a <- beta*2", "b <- a + gamma", "c <- b*a", "d <- c - b", "e <- d/(a+1)")
ode <- function(perm) {
  b <- compfit:::.build_model(mk(fn[perm]), backend = "r")
  b$compartmental_function
}
set.seed(11)
X <- c(990, 10); p <- c(beta = 0.3, gamma = 0.11)
ref <- ode(seq_along(fn))(X, p, 3)
chk("the reference build produces finite derivatives", all(is.finite(ref)))
same <- vapply(1:8, function(k) {
  identical(ode(sample(length(fn)))(X, p, 3), ref)
}, logical(1))
chk("8 random permutations give bitwise-identical derivatives", all(same))

th_section("empty / absent Functions columns are no-ops")
chk_equal("empty vector", .function_order(character(0))$order, integer(0))
chk_equal("all-blank column", .function_order(c(NA, "", "   "))$order, integer(0))

th_summary("function-order")
})
