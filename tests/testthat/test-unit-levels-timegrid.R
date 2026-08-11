test_that("test-unit-levels-timegrid", {
# ============================================================
# test-unit-levels-timegrid.R   (pure R; no Julia)
# Two behaviours:
#   * levels are detected by the '_' PREFIX (not the word "Level"), and a sheet
#     with NO level column defaults to a single level holding every compartment
#     (so N1 = total_pop and quadratic terms still normalise);
#   * .time_grid() warns when startpoint/endpoint is not a plain integer year
#     (a date or decimal is silently truncated to its leading number).
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R",
               "fitCompartmentalModel.R", "evaluate.R", "simulate.R"))

sc <- load_scenario(fixture_dir("SIR"), combined_file = "dataCombined.csv",
                    dummy_file = "dataDummy.csv", params_file = "modelParams.csv")

th_section("no level column -> single auto-level of all compartments")
mp0 <- sc$modelParams; mp0[["_Level1"]] <- NULL
chk("numberOfComps still counts via States", numberOfComps(mp0)$number_of_comps == 3)
chk("no '_' column detected", length(numberOfComps(mp0)$compartment_cols) == 0)
m0 <- chk_ok("quadratic model builds without a level column",
             build_compartmental_model(mp0, sc$dataCombined, solver = solver_control(backend = "r")))
chk("N1 spans every compartment",
    grepl("real N1 = X[1]+X[2]+X[3]", m0$model$stan_code, fixed = TRUE))

th_section("level column name is free (detected by '_' prefix, not \"Level\")")
mpA <- sc$modelParams; mpA[["_Age1"]] <- mpA[["_Level1"]]; mpA[["_Level1"]] <- NULL
chk("'_Age1' is recognised as a level column",
    identical(numberOfComps(mpA)$compartment_cols, "_Age1"))
chk_ok("model with _Age1 builds",
       build_compartmental_model(mpA, sc$dataCombined, solver = solver_control(backend = "r")))

th_section("no-level model simulates and conserves a closed population")
mps <- read_data_file(file.path(fixture_dir("SIR_sim"), "modelParams.csv"))
mps[["_Level1"]] <- NULL
dds <- read_data_file(file.path(fixture_dir("SIR_sim"), "dataDummy.csv"))
sim <- chk_ok("no-level SIR simulates",
              simulate_model(mps, data_dummy = dds, solver = solver_control(backend = "r")))
tot <- rowSums(sim$sir_out[c("X1", "X2", "X3")])
chk_equal("closed population conserved (N0 = 5000)", max(abs(tot - 5000)), 0, tol = 1e-2)

th_section("Mixing_<level> weights the transmission denominator (full pool)")
add_col <- function(mp, name, vals)
  { mp[[name]] <- c(as.character(vals), rep("", nrow(mp) - length(vals))); mp }
# Down-weight X2 in level 1's mixing pool: Nw1 = X1 + 0.5*X2 + X3.
mpM <- add_col(sc$modelParams, "Mixing_Level1", c("1", "0.5", "1"))
mM  <- chk_ok("model with Mixing_Level1 builds",
              build_compartmental_model(mpM, sc$dataCombined, solver = solver_control(backend = "r")))
chk("raw N1 (true head count) is left intact",
    grepl("real N1 = X[1]+X[2]+X[3];", mM$model$stan_code, fixed = TRUE))
chk("Stan defines the weighted pool Nw1",
    grepl("real Nw1 = X[1]+(0.5)*X[2]+X[3];", mM$model$stan_code, fixed = TRUE))
chk("Stan transmission divides by Nw1",
    grepl("/Nw1", mM$model$stan_code, fixed = TRUE))
chk("Julia defines Nw1 and divides by it",
    grepl("Nw1 = X[1]+(0.5)*X[2]+X[3]", mM$model$julia_code, fixed = TRUE) &&
    grepl("/Nw1", mM$model$julia_code, fixed = TRUE))
chk("R model body defines Nw1 and divides by it", {
  b <- paste(deparse(body(mM$model$compartmental_function)), collapse = "\n")
  grepl("Nw1 = X[1] + (0.5) * X[2] + X[3]", b, fixed = TRUE) && grepl("/Nw1", b, fixed = TRUE) })

th_section("no Mixing column -> byte-identical N-denominators (no Nw)")
m0m <- build_compartmental_model(sc$modelParams, sc$dataCombined, solver = solver_control(backend = "r"))
chk("no Nw pool is generated", !grepl("Nw", m0m$model$stan_code, fixed = TRUE))
chk("transmission still divides by the raw N1", grepl("/N1", m0m$model$stan_code, fixed = TRUE))

th_section("Mixing can pull a compartment from another level (cross-level pool)")
mp2 <- sc$modelParams
mp2[["_Level1"]] <- c("1", "2", rep("", nrow(mp2) - 2))   # level 1 = {X1, X2}
mp2 <- add_col(mp2, "_Level2", c("3"))                     # level 2 = {X3}
mp2 <- add_col(mp2, "Mixing_Level1", c("1", "1", "0.3"))   # X3 (level 2) enters level-1 pool
m2  <- chk_ok("two-level model with cross-level Mixing builds",
              build_compartmental_model(mp2, sc$dataCombined, solver = solver_control(backend = "r")))
chk("Nw1 includes the out-of-level compartment X[3]",
    grepl("real Nw1 = X[1]+X[2]+(0.3)*X[3];", m2$model$stan_code, fixed = TRUE))

th_section("Pool_<level>: the pool as a function of the level head counts")
mpP <- sc$modelParams
mpP[["_Level1"]] <- c("1", "2", rep("", nrow(mpP) - 2))   # level 1 = {X1, X2}
mpP <- add_col(mpP, "_Level2", c("3"))                    # level 2 = {X3}
mpP <- add_col(mpP, "Pool_Level1", c("N1 + 0.3*N2"))      # saturating cross-level pool
mP  <- chk_ok("Pool_<level> model builds",
              suppressMessages(build_compartmental_model(mpP, sc$dataCombined, solver = solver_control(backend = "r"))))
chk("Stan floors the pool with fmax (Stan's binary max)",
    grepl("real Nw1 = fmax((N1+0.3*N2), 1e-8);", mP$model$stan_code, fixed = TRUE))
chk("Julia floors the pool with max",
    grepl("Nw1 = max((N1+0.3*N2), 1e-8)", mP$model$julia_code, fixed = TRUE))
chk("R floors the pool with max", {
  b <- paste(deparse(body(mP$model$compartmental_function)), collapse = "\n")
  grepl("Nw1 = max((N1 + 0.3 * N2), 1e-08)", b, fixed = TRUE) })
chk("transmission divides by the pool Nw1", grepl("/Nw1", mP$model$stan_code, fixed = TRUE))

th_section("named level head-count aliases N_<level> (all backends)")
mpA <- sc$modelParams
mpA[["_Level1"]] <- NULL
mpA <- add_col(mpA, "_LA", c("1", "2"))   # level 1 -> alias N_LA
mpA <- add_col(mpA, "_HA", c("3"))        # level 2 -> alias N_HA
mpA$Functions[2] <- "share <- N_LA/(N_LA + N_HA)"   # reference levels BY NAME
mA <- chk_ok("model referencing N_<level> by name builds",
             build_compartmental_model(mpA, sc$dataCombined, solver = solver_control(backend = "r")))
chk("Stan emits the aliases and the by-name reference",
    grepl("real N_LA = N1;", mA$model$stan_code, fixed = TRUE) &&
    grepl("real N_HA = N2;", mA$model$stan_code, fixed = TRUE) &&
    grepl("share=N_LA/(N_LA+N_HA)", mA$model$stan_code, fixed = TRUE))
chk("Julia emits the aliases", grepl("N_LA = N1", mA$model$julia_code, fixed = TRUE))
chk("R body emits the aliases", {
  b <- paste(deparse(body(mA$model$compartmental_function)), collapse = "\n")
  grepl("N_LA = N1", b, fixed = TRUE) && grepl("N_HA = N2", b, fixed = TRUE) })

th_section("state names in Functions / Constant are rewritten to X[k] (all backends)")
mpN <- sc$modelParams
mpN$States   <- sub("X1", "S", mpN$States); mpN$States <- sub("X2", "I", mpN$States)
mpN$States   <- sub("X3", "R", mpN$States)
mpN$Functions[2] <- "prev <- I/(S+I+R)"                       # references states by name
mpN$Constant     <- c("0", "0.02*R", "0", rep("", nrow(mpN) - 3))  # inflow into I ~ R
mN <- chk_ok("model with state names in expressions builds",
             build_compartmental_model(mpN, sc$dataCombined, solver = solver_control(backend = "r")))
chk("Stan rewrites the Functions state refs",
    grepl("real prev=X[2]/(X[1]+X[2]+X[3]);", mN$model$stan_code, fixed = TRUE))
chk("Stan rewrites the Constant state ref (0.02*R -> 0.02*X[3])",
    grepl("0.02*X[3]", mN$model$stan_code, fixed = TRUE))
chk("Julia rewrites the same",
    grepl("prev=X[2]/(X[1]+X[2]+X[3])", mN$model$julia_code, fixed = TRUE))
chk("R model body rewrites the same", {
  b <- paste(deparse(body(mN$model$compartmental_function)), collapse = "\n")
  grepl("prev = X[2]/(X[1] + X[2] + X[3])", b, fixed = TRUE) })

th_section(".time_grid warns on a non-integer year")
warns <- function(mp) {
  w <- character(0)
  withCallingHandlers(compfit:::.time_grid(mp),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  w
}
set_start <- function(v) { mp <- sc$modelParams
  mp$Others[mp$Others == "startpoint=2010"] <- paste0("startpoint=", v); mp }
chk("date startpoint warns",    any(grepl("integer year", warns(set_start("2010-06-01")))))
chk("decimal startpoint warns", any(grepl("integer year", warns(set_start("2010.5")))))
chk("clean integer startpoint is silent", length(warns(sc$modelParams)) == 0)

th_summary("levels-timegrid")
})
