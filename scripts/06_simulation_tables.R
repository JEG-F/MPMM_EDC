###############################################################################
# 06_simulation_tables.R
#
#
# Tables 4, 6, 8 and 9 of the paper: proportion of misclassified pairs,
# Equation (12),
#     sum_r (Type 1 errors_r + Type 2 errors_r) / (1000 * choose(192, 2)),
# for every model, penalty and sample size N of Settings A-D.
#
# The input files hold the total numbers of Type 1 and Type 2 errors over
# the 1000 replications, as written by scripts/02-05:
#   results/simulation/setting_{A,B,C,D}_errors_seed171266.csv
#
# Usage (from the repository root):
#   Rscript scripts/06_simulation_tables.R                 # stored results (paper)
#   Rscript scripts/06_simulation_tables.R --input=results/simulation/reruns
#                                                          # results of a full re-run
#  Output: results/simulation/tables/table_setting_{A,B,C,D}.tex (LaTeX, as in
#         the paper) and results/simulation/tables/proportions.csv
#         (or, with --input, the same files inside <input>/tables/).
###############################################################################

root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else normalizePath(".")
})
if (!file.exists(file.path(root, "R", "edc_partition.R")))
  stop("Run this script from the repository root: Rscript scripts/06_simulation_tables.R")

args <- commandArgs(trailingOnly = TRUE)
in_arg <- sub("^--input=", "", grep("^--input=", args, value = TRUE))
in_dir <- if (length(in_arg)) file.path(root, in_arg) else file.path(root, "results", "simulation")
out_dir <- file.path(in_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed    <- 171266L
n_rep   <- 1000L
n_pairs <- choose(192, 2)           # pairs of elements of {1,...,6} x {0,1}^5
sizes   <- seq(1000L, 9000L, by = 1000L)
pens    <- c(w1 = "ln(n)", w2 = "(ln(ln(n)))^2", w3 = "ln(ln(n))")

# Models reported in each table, with the equations that define them
tables <- list(
  A = list(table = 4, label = "tab:simulation1",
           models = c("Model 1 (5 parts)", "Model 2 (10 parts)", "Model 3 (15 parts)"),
           head   = c("\\textbf{\\eqref{SAM1}}", "\\textbf{\\eqref{SAM2}}", "\\textbf{\\eqref{SAM3}}"),
           eqs    = c("(13)", "(14)", "(15)")),
  B = list(table = 6, label = "tab:simulation2",
           models = c("Model 2 (10 parts)", "Model 3 (15 parts)"),
           head   = c("\\textbf{\\eqref{SBM1}}", "\\textbf{\\eqref{SBM2}}"),
           eqs    = c("(16)", "(17)")),
  C = list(table = 8, label = "tab:simulation3",
           models = c("Model 2 (10 parts)", "Model 3 (15 parts)"),
           head   = c("\\textbf{\\eqref{SCM1a}}-\\textbf{\\eqref{SCM1c}}",
                      "\\textbf{\\eqref{SCM2a}}-\\textbf{\\eqref{SCM2d}}"),
           eqs    = c("(18)-(20)", "(21)-(24)")),
  D = list(table = 9, label = "tab:simulation4",
           models = c("Model 1 (10 parts)", "Model 2 (10 parts)", "Model 3 (10 parts)"),
           head   = c("\\textbf{\\eqref{SDM1}}", "\\textbf{\\eqref{SDM2}}", "\\textbf{\\eqref{SDM3}}"),
           eqs    = c("(25)", "(26)", "(27)"))
)

all_props <- NULL

for (s in names(tables)) {
  tb <- tables[[s]]
  f  <- file.path(in_dir, sprintf("setting_%s_errors_seed%d.csv", s, seed))
  if (!file.exists(f)) {
    cat("Setting", s, ": file not found, skipped:", f, "\n")
    next
  }
  res <- read.csv(f, stringsAsFactors = FALSE)
  res$prop <- (res$Type1 + res$Type2) / (n_rep * n_pairs)

  # Matrix of proportions: rows N, columns model x penalty (as in the paper)
  props <- sapply(seq_along(tb$models), function(i) sapply(pens, function(pn)
    sapply(sizes, function(N) {
      v <- res$prop[res$Model == tb$models[i] & res$Penalty == pn & res$N == N]
      if (length(v) == 1) v else NA
    })), simplify = "array")
  props <- matrix(props, nrow = length(sizes))   # columns: model 1 (w1, w2, w3), model 2 (...), ...
  if (anyNA(props)) stop("Setting ", s, ": incomplete result file ", f)

  all_props <- rbind(all_props, data.frame(
    setting  = s,
    table    = tb$table,
    model    = rep(tb$models, each = 3 * length(sizes)),
    equation = rep(tb$eqs, each = 3 * length(sizes)),
    penalty  = rep(rep(names(pens), each = length(sizes)), length(tb$models)),
    N        = rep(sizes, 3 * length(tb$models)),
    proportion = as.numeric(sprintf("%.4f", as.vector(props))),
    stringsAsFactors = FALSE))

# Plain-text version
  cat(sprintf("\nTable %d (Setting %s): proportion of errors, Equation (12)\n", tb$table, s))
  hdr <- paste(sprintf("%-26s", paste0("Eq. ", tb$eqs)), collapse = "")
  cat(sprintf("%6s  %s\n", "", hdr))
  cat(sprintf("%6s  %s\n", "N", paste(rep(sprintf("%-8s %-8s %-8s ", "w1", "w2", "w3"),
                                           length(tb$models)), collapse = " ")))
  for (r in seq_along(sizes))
    cat(sprintf("%6d  %s\n", sizes[r], paste(sprintf("%.4f  ", props[r, ]), collapse = " ")))

  # LaTeX version (same layout as in the paper)
  nc  <- 3 * length(tb$models)
  out <- c(sprintf("\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}l%s@{}}", strrep("c", nc)),
           "\\toprule",
           paste0("\\textbf{Model} & ",
                  paste(sprintf("\\multicolumn{3}{c}{%s}", tb$head), collapse = " & "), " \\\\"),
           paste(sprintf("\\cmidrule{%d-%d}", 3 * seq_along(tb$models) - 1,
                         3 * seq_along(tb$models) + 1), collapse = " "),
           paste0("\\textbf{$w$} & ",
                  paste(rep("\\textbf{$w_1$} & \\textbf{$w_2$} & \\textbf{$w_3$}",
                            length(tb$models)), collapse = " & "), " \\\\"),
           "\\midrule",
           paste0("$N$", strrep(" &", nc), " \\\\"),
           sprintf("%d & %s \\\\", sizes,
                   apply(props, 1, function(x) paste(sprintf("%.4f", x), collapse = " & "))),
           "\\bottomrule",
           "\\end{tabular*}")
  writeLines(out, file.path(out_dir, sprintf("table_setting_%s.tex", s)))
}
write.csv(all_props, file.path(out_dir, "proportions.csv"), row.names = FALSE)
cat("\nTables written to:", out_dir, "\n")
