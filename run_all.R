###############################################################################
# run_all.R -- reproduce the results of the paper
# Rscript run_all.R
#   Exchange-rate application (Tables 10-16, Figures 3-5) and simulation
#       tables (Tables 4, 6, 8, 9) from the stored error totals, then every
#       number is compared with the value printed in the paper.
#
#   Rscript run_all.R --simulations --cores=8
#       In addition, it re-runs the four Monte Carlo studies (Settings A-D) on
#       8 cores, compares every cell with the stored totals and rebuilds the
#       tables from the re-run (results/simulation/reruns/).  This may takes hours,
#       see README.md.
###############################################################################

root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) normalizePath(dirname(f)) else normalizePath(".")
})
if (!file.exists(file.path(root, "R", "edc_partition.R")))
  stop("Run this script from the repository root: Rscript run_all.R")

args      <- commandArgs(trailingOnly = TRUE)
run_sims  <- "--simulations" %in% args
cores_arg <- grep("^--cores=[0-9]+$", args, value = TRUE)
rscript   <- file.path(R.home("bin"), "Rscript")

run <- function(script, ...) {
  cat("\n", strrep("#", 78), "\n# Rscript ", script, " ", paste(...), "\n",
      strrep("#", 78), "\n", sep = "")
  status <- system2(rscript, c(shQuote(file.path(root, script)), ...))
  if (status != 0) stop(script, " failed (exit status ", status, ")")
}

run("scripts/01_application_exchange_rates.R")

if (run_sims) {
  for (s in c("02_simulation_setting_A.R", "03_simulation_setting_B.R",
              "04_simulation_setting_C.R", "05_simulation_setting_D.R"))
    run(file.path("scripts", s), cores_arg)
  run("scripts/06_simulation_tables.R", "--input=results/simulation/reruns")
}

run("scripts/06_simulation_tables.R")

invisible(requireNamespace("ape", quietly = TRUE))   # record its version below
writeLines(capture.output(sessionInfo()), file.path(root, "results", "sessionInfo.txt"))

run("tests/check_against_paper.R")
