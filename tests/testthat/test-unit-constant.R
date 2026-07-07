test_that("test-unit-constant", {
# ============================================================
# test-unit-constant.R   (pure R; no Julia)
# The optional single `Constant` column: row i is a zeroth-order term added
# directly to dX[i] (no state multiplier) -- e.g. a constant inflow / birth
# rate. Verifies it (a) lands in the per-compartment derivative on all three
# code backends (R body, Julia, Stan), (b) accepts a number and a parameter
# name, (c) an absent/blank column adds nothing (silently), and (d) the
# validator recognises the column and still flags an unknown symbol in a cell.
# ============================================================

th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R", "validate.R"))

load_ex <- function(name)
  load_scenario(fixture_dir(name), combined_file = "dataCombined.csv",
                dummy_file = "dataDummy.csv", params_file = "modelParams.csv")

dX1 <- function(v) grep("dX\\[?1\\]?\\s*=", v, value = TRUE)[1]
body_lines <- function(mod)
  strsplit(paste(deparse(body(mod$compartmental_function)), collapse = "\n"), "\n")[[1]]

sc <- load_ex("SIR")
mp <- sc$modelParams
mp$Constant    <- rep("", nrow(mp))
mp$Constant[1] <- "1234"      # numeric inflow into compartment 1
mp$Constant[3] <- "gamma"     # parameter-valued inflow into compartment 3
mod <- suppressWarnings(compfit:::.build_model(mp, backend = "r"))

rlines <- body_lines(mod)
jlines <- strsplit(mod$julia_code, "\n")[[1]]
slines <- strsplit(mod$stan_code,  "\n")[[1]]

th_section("Constant term appears in the derivative on every backend")
chk("R    dX1 carries the numeric constant", grepl("1234", dX1(rlines)))
chk("Julia dX1 carries the numeric constant", grepl("1234", dX1(jlines)))
chk("Stan  dX1 carries the numeric constant", grepl("1234", dX1(slines)))
chk("parameter-valued constant lands in dX3 (R)",
    grepl("gamma", grep("dX3", rlines, value = TRUE)[1]))

th_section("absent Constant column -> no constant term (silent)")
mod0 <- suppressWarnings(compfit:::.build_model(sc$modelParams, backend = "r"))
r0   <- body_lines(mod0)
chk("sheet without a Constant column still builds", !is.null(mod0$compartmental_function))
chk("no stray constant when the column is absent", !grepl("1234", paste(r0, collapse = " ")))

th_section("blank Constant cell -> 0 (no term for that compartment)")
mpb <- mp; mpb$Constant[1] <- ""            # blank the compartment-1 constant
modb <- suppressWarnings(compfit:::.build_model(mpb, backend = "r"))
chk("blanked constant cell drops the term", !grepl("1234", dX1(body_lines(modb))))

th_section("validator recognises Constant and still checks its cells")
w  <- character(0)
ok <- withCallingHandlers(validate_modelParams(mp),
        warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
chk("valid Constant column -> no unrecognised-column warning", !any(grepl("unrecognised", w)))
mpbad <- mp; mpbad$Constant[1] <- "notaparam"
chk("unknown symbol in a Constant cell is flagged",
    isTRUE(tryCatch({ validate_modelParams(mpbad); FALSE },
                    error = function(e) grepl("Constant", conditionMessage(e)))))

th_summary("constant")
})
