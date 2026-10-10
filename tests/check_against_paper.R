###############################################################################
# check_against_paper.R
#
# To simplify the checking process, this script compares the outputs of the 
# scripts with the numbers printed in the paper.
# The expected values in tests/expected/ were extracted from the LaTeX source
# of the manuscript (version of 2026-10-08):
#   paper_simulation_tables.csv  Tables 4, 6, 8 and 9 (270 proportions)
#   paper_table10.csv            Table 10  (LML, EDC, |L|)
#   paper_partitions.csv         Tables 11-12 (w1 and w2) and 14-15 (w3):
#                                cardinality, P(0|L), P(1|L) and the members
#                                of every part, in the printed order
#   paper_table13.csv            Table 13 (states per currency)
#   paper_table16.csv            Table 16 (w1 parts split into w3 parts)
#
# Run scripts/01_application_exchange_rates.R and
#  scripts/06_simulation_tables.R first (run_all.R does both).
#
#Usage (from the repository root): Rscript tests/check_against_paper.R
# Exit status 1 if any check fails.
###############################################################################

root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else normalizePath(".")
})
exp_dir <- file.path(root, "tests", "expected")
app_dir <- file.path(root, "results", "application")

all_ok <- TRUE
report <- function(what, ok, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", what,
              if (nzchar(detail)) paste0("  (", detail, ")") else ""))
  all_ok <<- all_ok && isTRUE(ok)
}
need <- function(f) {
  if (!file.exists(f)) {
    report(basename(f), FALSE, "file not found; run run_all.R first")
    return(FALSE)
  }
  TRUE
}

cat("Checking the outputs against the values printed in the paper\n\n")

# ---- Tables 4, 6, 8, 9 ----------------------------------------------------

f <- file.path(root, "results", "simulation", "tables", "proportions.csv")
if (need(f)) {
  got <- read.csv(f, stringsAsFactors = FALSE)
  exp <- read.csv(file.path(exp_dir, "paper_simulation_tables.csv"), stringsAsFactors = FALSE)
  key <- function(d) paste(d$setting, d$model, d$penalty, d$N)
  m <- match(key(exp), key(got))
  for (s in unique(exp$setting)) {
    i <- which(exp$setting == s)
    same <- !is.na(m[i]) & abs(got$proportion[m[i]] - exp$proportion[i]) < 1e-9
    same[is.na(same)] <- FALSE
    report(sprintf("Table %d (Setting %s): %d proportions", exp$table[i[1]], s, length(i)),
           all(same), if (!all(same)) sprintf("%d differ", sum(!same)) else "")
  }
}

# --- Table 10 ---------------------------------------------------------------

f <- file.path(app_dir, "summary_LML_EDC.csv")
if (need(f)) {
  got <- read.csv(f, stringsAsFactors = FALSE)
  exp <- read.csv(file.path(exp_dir, "paper_table10.csv"), stringsAsFactors = FALSE)
  ok <- identical(got$penalty, exp$penalty) &&
    all(abs(round(got$LML, 2) - exp$LML) < 1e-9) &&
    all(abs(round(got$EDC, 2) - exp$EDC) < 1e-9) &&
    all(got$parts == exp$parts)
  report("Table 10: LML, EDC and |L| for the three penalties", ok)
}

# --- Tables 11-12, 14-15, 13 and 16 --------------------------------------------

f <- file.path(app_dir, "partitions.csv")
if (need(f)) {
  got <- read.csv(f, colClasses = c(state = "character"), stringsAsFactors = FALSE)
  exp <- read.csv(file.path(exp_dir, "paper_partitions.csv"),
                  colClasses = c(state = "character"), stringsAsFactors = FALSE)

  tables <- c("ln(n)" = "Tables 11-12 (w1 = ln(n))",
              "(ln(ln(n)))^2" = "Tables 11-12 (w2 = [ln(ln(n))]^2)",
              "ln(ln(n))" = "Tables 14-15 (w3 = ln(ln(n)))")
  for (pen in names(tables)) {
    g <- got[got$penalty == pen, ]
    e <- exp[exp$penalty == pen, ]
    parts <- sort(unique(e$part))
    ok <- identical(sort(unique(g$part)), parts)
    if (ok) for (v in parts) {
      gv <- g[g$part == v, ]
      ev <- e[e$part == v, ]
      ok <- ok && nrow(gv) == nrow(ev) && gv$n_elements[1] == ev$n_elements[1] &&
        abs(round(gv$P0[1], 4) - ev$P0[1]) < 1e-9 &&
        abs(round(gv$P1[1], 4) - ev$P1[1]) < 1e-9 &&
        identical(paste(gv$currency, gv$state), paste(ev$currency[order(ev$order)], ev$state[order(ev$order)]))
    }
    report(sprintf("%s: %d parts, all members", tables[[pen]], length(parts)), ok)
  }

  # Table 13
  e13 <- read.csv(file.path(exp_dir, "paper_table13.csv"), stringsAsFactors = FALSE)
  g1  <- got[got$penalty == "ln(n)", ]
  cur <- setdiff(names(e13), c("part", "n_elements"))
  g13 <- t(sapply(e13$part, function(v) table(factor(g1$currency[g1$part == v], levels = cur))))
  report("Table 13: states per currency in each part",
         all(g13 == as.matrix(e13[, cur])) &&
           all(rowSums(g13) == e13$n_elements))

  # Table 16
  e16 <- read.csv(file.path(exp_dir, "paper_table16.csv"), stringsAsFactors = FALSE)
  g3  <- got[got$penalty == "ln(ln(n))", ]
  key1 <- setNames(g1$part, paste(g1$currency, g1$state))
  ok <- TRUE
  for (r in seq_len(nrow(e16))) {
    sub3 <- g3[g3$part == e16$part_w3[r], ]
    ok <- ok && nrow(sub3) == e16$n_elements_w3[r] &&
      all(key1[paste(sub3$currency, sub3$state)] == e16$part_w1[r]) &&
      abs(round(sub3$P0[1], 4) - e16$P0_w3[r]) < 1e-9
  }
  report("Table 16: decomposition of the w1 parts into the w3 parts", ok)
}

cat(if (all_ok) "\nAll checks passed.\n" else "\nSome checks FAILED.\n")
quit(status = if (all_ok) 0 else 1, save = "no")
