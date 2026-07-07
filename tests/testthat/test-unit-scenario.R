test_that("test-unit-scenario", {
# ============================================================
# test-unit-scenario.R   (pure R; no Julia)
# next_plot_path numbering, read_data_file text_cols coercion, save_grid no-op,
# read_data_file unsupported-extension error.
# ============================================================
th_load_pure(c("scenario.R"))

th_section("next_plot_path: auto-incrementing numbering")
# tempfile() gives a fresh, unique path on every call -- independent of the RNG
# seed and safe to re-run in the same session (a runif-based name could collide
# with a previous run's leftover Plot_*.pdf and break the numbering check).
d <- tempfile("npp_")
dir.create(d)
chk("first save -> _1", basename(next_plot_path(file.path(d, "Plot.pdf"))) == "Plot_1.pdf")
invisible(file.create(file.path(d, "Plot_3.pdf")))
chk("max+1 after Plot_3 -> _4", basename(next_plot_path(file.path(d, "Plot.pdf"))) == "Plot_4.pdf")
invisible(file.create(file.path(d, "Plot.pdf")))            # legacy un-numbered
chk("legacy Plot.pdf ignored -> still _4", basename(next_plot_path(file.path(d, "Plot.pdf"))) == "Plot_4.pdf")
chk("stem with regex chars handled", {
  file.create(file.path(d, "a.b_2.pdf"))
  basename(next_plot_path(file.path(d, "a.b.pdf"))) == "a.b_3.pdf"
})

th_section("read_data_file: text_cols preserves stored-as-text numbers")
csv <- file.path(d, "tc.csv")
writeLines(c("a,b", "0,1", "2,3"), csv)
num <- read_data_file(csv)                 # normal typing
chk("default read is numeric", is.numeric(num$a) && num$a[1] == 0)
txt <- read_data_file(csv, text_cols = TRUE)
chk("text_cols read is character", is.character(txt$a))
chk("text '0' survives as character", txt$a[1] == "0")
chk_equal("as.numeric round-trips", as.numeric(txt$a), c(0, 2))

th_section("read_data_file: unsupported extension errors")
invisible(file.create(file.path(d, "x.weird")))
chk_error("unknown extension", read_data_file(file.path(d, "x.weird")))

th_section("read_data_file: read failures name the file")
# a file that is not a real .xlsx -> error must mention the path (which file failed)
bx <- file.path(d, "broken.xlsx"); writeLines("not an excel file", bx)
msg <- tryCatch({ read_data_file(bx); "" }, error = function(e) conditionMessage(e))
chk("broken read names the file", grepl("broken.xlsx", msg, fixed = TRUE))
chk("broken read says 'Could not read'", grepl("Could not read", msg, fixed = TRUE))
# missing file -> clear error naming the path
miss <- file.path(d, "nope.csv")
chk_error("missing file errors", read_data_file(miss))
chk("missing-file error names the path",
    grepl("nope.csv", tryCatch({ read_data_file(miss); "" }, error = function(e) conditionMessage(e)), fixed = TRUE))

th_section("save_grid: no-grid branch is a safe no-op")
chk("returns NULL when $grid is NULL", is.null(save_grid(list(grid = NULL), file.path(d, "P.pdf"))))

th_summary("scenario")
})
