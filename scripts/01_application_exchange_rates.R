###############################################################################
# 01_application_exchange_rates.R
#
# Section 4 of the paper: daily direction of six exchange rates against the
# US dollar (BRL, EUR, CAD, GBP, MXN, CNY), 2006-03-24 to 2026-02-02.
# Binary alphabet A = {0,1}, Markov order o = 5, p = 6; Algorithm 1 with
# gamma = identity, alpha = 2 and the penalties w1 = ln(n),
# w2 = [ln(ln(n))]^2 and w3 = ln(ln(n)).
# (Original script: run_application.R.)
#
# Reproduces
#   Table 10      LML, EDC and |L| for the three penalties
#   Tables 11-12  estimated partition with w1 (identical with w2)
#   Table 13      number of states of each currency in each part (w1, w2)
#   Tables 14-15  estimated partition with w3
#   Table 16      decomposition of the w1/w2 parts into the w3 parts
#   Figures 3-5   circular dendrograms of d for w1, w2 and w3
#
# Input : data/dados_yahoo.RDATA
# Output: data/binary_series.csv                  binary series (5182 x 6)
#         results/application/console_output.txt  copy of the console output
#         results/application/summary_LML_EDC.csv Table 10
#         results/application/partitions.csv      Tables 11-12 and 14-15
#         results/application/tables/*.tex        LaTeX of Tables 10-16
#         results/application/figures/*.pdf       dendrograms; Figures 3-5 are
#             dendrograms_application_circular-1.pdf, -2.pdf and -3.pdf
#
# Usage (from the repository root):
#   Rscript scripts/01_application_exchange_rates.R
# Running time: a few seconds.  Requires the package 'ape'.
###############################################################################


root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else normalizePath(".")
})
if (!file.exists(file.path(root, "R", "edc_partition.R")))
  stop("Run this script from the repository root: Rscript scripts/01_application_exchange_rates.R")

source(file.path(root, "R", "edc_partition.R"))
source(file.path(root, "R", "edc_dendrogram.R"))
library(ape)


app_dir <- file.path(root, "results", "application")
tab_dir <- file.path(app_dir, "tables")
fig_dir <- file.path(app_dir, "figures")
for (d in c(tab_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Keep a copy of everything printed below
log_con <- file(file.path(app_dir, "console_output.txt"), open = "wt")
sink(log_con, split = TRUE)

# --- Data preparation --------------------------------------------------------
# dados_yahoo hold the daily quotes downloaded from Yahoo Finance; its column
# 'sinal' is the binary series of the paper (1 = increase or maintenance,
# 0 = decrease).See data/README.md.

load(file.path(root, "data", "dados_yahoo.RDATA"))

# Short labels for the currencies (this is also the column order j = 1..6
# used by the algorithm; the estimated partition does not depend on it)
currency_map <- c(
  "BRLUSD=X" = "BRL",
  "EURUSD=X" = "EUR",
  "CADUSD=X" = "CAD",
  "GBPUSD=X" = "GBP",
  "MXNUSD=X" = "MXN",
  "CNYUSD=X" = "CNY"
)

# Pivot to wide format: one column per currency, rows = trading days
dates <- sort(unique(dados_yahoo$date))

wide <- data.frame(row.names = seq_along(dates))
for (sym in names(currency_map)) {
  sub <- dados_yahoo[dados_yahoo$symbol == sym, ]
  sub <- sub[order(sub$date), ]
  stopifnot(nrow(sub) == length(dates), all(sub$date == dates))
  wide[[currency_map[sym]]] <- as.character(as.integer(sub$sinal))
}

binary_file <- file.path(root, "data", "binary_series.csv")
write.csv(wide, binary_file, row.names = FALSE)
cat("Binary series written to: data/binary_series.csv\n")
cat("  Rows (trading days):", nrow(wide), "\n")
cat("  Columns (currencies):", ncol(wide), "\n\n")

# --- EDC analysis ------------------------------------------------------------

order_o <- 5L

data <- read.csv(binary_file, header = TRUE, colClasses = "character")
counts_info <- compute_counts(data, order_o)
cn <- colnames(data)

cat("=== EDC Application: Exchange Rate Data ===\n\n")
cat("Currencies (p):  ", counts_info$p, "\n")
cat("Sample size (N): ", counts_info$n[1], "\n")
cat("Alphabet:         ", paste(counts_info$alphabet, collapse = ", "), "\n")
cat("Markov order (o): ", order_o, "\n")
cat("States |A^o|:     ", length(counts_info$states), "\n")
cat("Active (j,s):     ", sum(counts_info$state_counts > 0), "\n\n")

# --- Penalty functions (Table 1) ---------------------------------------------

penalties <- list(
  list(name = "BIC: ln(n)",         short = "ln(n)",
       tex  = "$\\ln(\\cdot)$",
       w_fn = function(n) log(n)),
  list(name = "EDC: (ln(ln(n)))^2", short = "(ln(ln(n)))^2",
       tex  = "$[\\ln(\\ln(\\cdot))]^2$",
       w_fn = function(n) log(log(n))^2),
  list(name = "EDC: ln(ln(n))",     short = "ln(ln(n))",
       tex  = "$\\ln(\\ln(\\cdot))$",
       w_fn = function(n) log(log(n)))
)

# --- Run each penalty --------------------------------------------------------

fits <- list()
summary_rows <- list()

for (k in seq_along(penalties)) {
  pen <- penalties[[k]]
  cat(strrep("=", 70), "\n")
  cat("  Penalty:", pen$name, "\n")
  cat(strrep("=", 70), "\n\n")

  result <- edc_partition(counts_info, gamma_fn = identity,
                          w_fn = pen$w_fn, alpha = 2)
  fits[[k]] <- result

  cat("  Parts:", length(result$partition),
      " | Merges:", length(result$merge_history), "\n\n")

  for (v in seq_along(result$partition)) {
    members <- result$partition[[v]]
    cat("  Part", v, "(", nrow(members), "elements )\n")
    cat("    P(a|L):", paste(paste0(counts_info$alphabet, "=",
                                     round(result$transition_probs[v, ], 4)),
                              collapse = "  "), "\n")
    cat("    Members:")
    for (i in seq_len(nrow(members))) {
      # Use currency name instead of column number
      cat(" (", cn[members$sample[i]], ",", members$state[i], ")")
    }
    cat("\n")
  }

  # LML (Equation (6)) and EDC (Equation (7)) of the estimated partition
  summary_rows[[k]] <- data.frame(
    penalty = pen$short,
    LML     = result$lml_value,
    EDC     = result$edc_value,
    parts   = length(result$partition),
    stringsAsFactors = FALSE
  )

  cat("\n")
}

summary_tab <- do.call(rbind, summary_rows)

# ----- Summary table (Table 10) ------------------------------------------------

cat("\n")
cat(strrep("=", 60), "\n")
cat("  Summary (Table 10)\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  %-18s %12s %12s %6s\n", "w(n)", "LML", "EDC", "Parts"))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(summary_tab))) {
  cat(sprintf("  %-18s %12.2f %12.2f %6d\n", summary_tab$penalty[i],
              summary_tab$LML[i], summary_tab$EDC[i], summary_tab$parts[i]))
}
cat(strrep("-", 60), "\n\n")
write.csv(summary_tab, file.path(app_dir, "summary_LML_EDC.csv"), row.names = FALSE)

# ------ Partitions in long format (Tables 11-12 and 14-15) ----------------------

partitions <- do.call(rbind, lapply(seq_along(penalties), function(k) {
  res <- fits[[k]]
  do.call(rbind, lapply(seq_along(res$partition), function(v) {
    m <- res$partition[[v]]
    data.frame(penalty    = penalties[[k]]$short,
               part       = v,
               n_elements = nrow(m),
               P0         = unname(res$transition_probs[v, "0"]),
               P1         = unname(res$transition_probs[v, "1"]),
               currency   = cn[m$sample],
               state      = m$state,
               stringsAsFactors = FALSE)
  }))
}))
write.csv(partitions, file.path(app_dir, "partitions.csv"), row.names = FALSE)

# Membership label of every (currency, state) pair in a fitted partition
part_labels <- function(res) {
  lab <- integer(0)
  for (v in seq_along(res$partition)) {
    m <- res$partition[[v]]
    lab[paste0("(", cn[m$sample], ",", m$state, ")")] <- v
  }
  lab
}

same_w1_w2 <- identical(part_labels(fits[[1]]), part_labels(fits[[2]]))
cat("Partitions for w1 = ln(n) and w2 = [ln(ln(n))]^2 identical:", same_w1_w2, "\n\n")

# ----- LaTeX tables ------------------------------------------------------------

fmt4 <- function(x) formatC(round(x, 4), format = "f", digits = 4)
fmt_big <- function(x) gsub(",", "{,}", formatC(x, format = "f", digits = 2, big.mark = ","))

# Table 10
writeLines(c(
  "\\begin{tabular}{lrrr}",
  "\\toprule",
  "$w(\\cdot)$ & $\\mbox{LML}_{\\widehat{\\mathbb{L}}}$ & $\\mbox{EDC}_{\\widehat{\\mathbb{L}}}$ & $|\\widehat{\\mathbb{L}}|$  \\\\",
  "\\midrule",
  sprintf("%-23s & $%s$ & $%s$ & %d  \\\\", vapply(penalties, `[[`, "", "tex"),
          fmt_big(summary_tab$LML), fmt_big(summary_tab$EDC), summary_tab$parts),
  "\\bottomrule",
  "\\end{tabular}"), file.path(tab_dir, "summary_LML_EDC.tex"))

# Tables 11-12 and 14-15: one tabular per penalty (the paper splits each one
# over two tables for page layout); five members per line
partition_tex <- function(res) {
  out <- c("\\begin{tabular}{lllll}", "\\toprule",
           "Part($L$) & $N(L)$ & $P(0|L)$ & $P(1|L)$ & Members \\\\", "\\midrule")
  K <- length(res$partition)
  for (v in seq_len(K)) {
    m <- res$partition[[v]]
    mem <- sprintf("(%s,%s)", cn[m$sample], m$state)
    lines <- vapply(split(mem, ceiling(seq_along(mem) / 5)), paste, "", collapse = " ")
    out <- c(out,
             sprintf("$L_{%d}$ & %d & %s & %s &", v, nrow(m),
                     fmt4(res$transition_probs[v, 1]), fmt4(res$transition_probs[v, 2])),
             "  \\begin{minipage}[t]{0.55\\textwidth}\\raggedright\\scriptsize",
             paste0("  ", lines),
             paste0("  \\end{minipage} \\\\", if (v < K) "[2pt]" else ""))
  }
  c(out, "\\bottomrule", "\\end{tabular}")
}
writeLines(partition_tex(fits[[1]]), file.path(tab_dir, "partition_w1.tex"))
writeLines(partition_tex(fits[[2]]), file.path(tab_dir, "partition_w2.tex"))
writeLines(partition_tex(fits[[3]]), file.path(tab_dir, "partition_w3.tex"))

# Table 13: states of each currency in each part of the w1 (= w2) partition
by_currency <- t(vapply(fits[[1]]$partition, function(m)
  as.integer(table(factor(cn[m$sample], levels = cn))), integer(length(cn))))
colnames(by_currency) <- cn
by_currency <- cbind(by_currency, total = rowSums(by_currency))
rownames(by_currency) <- paste0("L", seq_len(nrow(by_currency)))
cat("States per currency in each part, w1 = ln(n) (Table 13):\n")
print(rbind(by_currency, Total = colSums(by_currency)))
cat("\n")
writeLines(c(
  "\\begin{tabular}{llllllll} \\toprule",
  paste0("Part  ($L$)&", paste(cn, collapse = "& "), "& $N(L)$\\\\ \\midrule"),
  sprintf("$L_{%d}$ &%s\\\\", seq_len(nrow(by_currency)),
          apply(by_currency, 1, paste, collapse = "& ")),
  "\\midrule",
  paste0("Total  by currency&", paste(colSums(by_currency), collapse = "& "), "\\\\"),
  "\\bottomrule",
  "\\end{tabular}"), file.path(tab_dir, "states_per_currency_w1.tex"))

# Table 16: each part of the w3 partition lies inside one part of the w1
# partition; rows ordered by increasing P(0|L)
lab1 <- part_labels(fits[[1]])
lab3 <- part_labels(fits[[3]])
parent <- lapply(seq_along(fits[[3]]$partition), function(v) unique(lab1[names(lab3)[lab3 == v]]))
if (all(lengths(parent) == 1)) {
  parent <- unlist(parent)
  p0_1 <- fits[[1]]$transition_probs[, 1]
  p0_3 <- fits[[3]]$transition_probs[, 1]
  size <- function(res, v) nrow(res$partition[[v]])
  out <- c("\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lccclcc@{}}", "\\toprule",
           "\\multicolumn{3}{c}{\\textbf{Estimated parts - $w_1$, $w_2$}} & & \\multicolumn{3}{c}{\\textbf{Estimated parts - $w_3$}} \\\\ ",
           "\\cmidrule{1-3} \\cmidrule{5-7} ",
           "\\textbf{Part ($L$)} & \\textbf{$N(L)$} & \\textbf{$P(0 \\vert L)$} & & \\textbf{Part ($L$)} & \\textbf{$N(L)$} & \\textbf{$P(0 \\vert L)$} \\\\ ",
           "\\midrule")
  cat("Decomposition of the w1 parts into the w3 parts (Table 16):\n")
  for (u in order(p0_1)) {
    kids <- which(parent == u)
    kids <- kids[order(p0_3[kids])]
    for (i in seq_along(kids)) {
      left <- if (i == 1) sprintf("$L_{%d}$ & %d & %s & $\\to$", u, size(fits[[1]], u), fmt4(p0_1[u]))
              else "      &    &        &      "
      end <- if (i == length(kids)) " \\\\ \\hline" else " \\\\"
      out <- c(out, sprintf("%s & $L_{%d}$ & %d & %s%s", left, kids[i],
                            size(fits[[3]], kids[i]), fmt4(p0_3[kids[i]]), end))
      cat(sprintf("  %-26s -> L%d (%2d, P0=%s)\n",
                  if (i == 1) sprintf("L%d (%2d, P0=%s)", u, size(fits[[1]], u), fmt4(p0_1[u])) else "",
                  kids[i], size(fits[[3]], kids[i]), fmt4(p0_3[kids[i]])))
    }
  }
  out[length(out)] <- sub(" \\\\hline$", "", out[length(out)])
  writeLines(c(out, "\\bottomrule", "\\end{tabular*}"),
             file.path(tab_dir, "partition_decomposition_w1_w3.tex"))
  cat("\n")
} else {
  cat("The w3 partition is not a refinement of the w1 partition; Table 16 not written.\n\n")
}

# ----- Dendrograms (Figures 3-5) -----------------------------------------------

trees <- lapply(penalties, function(pen)
  edc_full_tree(counts_info, cn, gamma_fn = identity, w_fn = pen$w_fn, alpha = 2))

plot_linear <- function(hc, pen) {
  par(mar = c(6, 4, 3, 1))
  n_clust <- sum(hc$height > 1) + 1
  plot(hc, hang = -1,
       main = paste0("Exchange Rate Application - w(n) = ", pen$short,
                     "  [", n_clust, " parts at d=1]"),
       xlab = "", ylab = "d",
       cex = 0.45, cex.main = 1.1)
  abline(h = 1, col = "red", lty = 2, lwd = 1.5)
  text(0.5, 1, "d = 1", col = "red", pos = 3, cex = 0.8)
}

plot_logscale <- function(hc, pen) {
  par(mar = c(6, 4, 3, 1))
  n_clust <- sum(hc$height > 1) + 1
  hc$height <- log(hc$height + 1)
  cut_transformed <- log(1 + 1)
  plot(hc, hang = -1,
       main = paste0("Exchange Rate Application - w(n) = ", pen$short,
                     "  [", n_clust, " parts at d=1]"),
       xlab = "", ylab = paste(expression(ln(d + 1)), " ", "scale"),
       cex = 0.45, cex.main = 1.1, axes = FALSE)
  max_h <- max(hc$height)
  d_ticks <- c(0, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
  d_ticks <- d_ticks[log(d_ticks + 1) <= max_h * 1.05]
  axis(2, at = log(d_ticks + 1), labels = d_ticks, las = 1, cex.axis = 0.7)
  box()
  abline(h = cut_transformed, col = "red", lty = 2, lwd = 1.5)
  text(0.5, cut_transformed, "d = 1", col = "red", pos = 3, cex = 0.8)
}

# Fan layout on the ln(d + 1) scale; the red dashed circle marks d = 1
plot_circular <- function(hc, pen) {
  par(mar = c(1, 1, 3, 1))
  n_clust <- sum(hc$height > 1) + 1
  hc$height <- log(hc$height + 1)
  phy <- as.phylo(hc)
  plot(phy, type = "fan",
       tip.color = "black", cex = 0.75,
       label.offset = 0.02,
       main = paste0("Exchange Rate Application - w(n) = ", pen$short,
                     "  [", n_clust, " parts at d=1]"))
  pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  n_tips  <- length(phy$tip.label)
  root_idx <- n_tips + 1
  root_x <- pp$xx[root_idx]
  root_y <- pp$yy[root_idx]
  tip_r  <- mean(sqrt((pp$xx[1:n_tips] - root_x)^2 +
                      (pp$yy[1:n_tips] - root_y)^2))
  max_h  <- max(hc$height)
  cut_r  <- tip_r * (max_h - log(1 + 1)) / max_h
  theta  <- seq(0, 2 * pi, length.out = 360)
  lines(root_x + cut_r * cos(theta), root_y + cut_r * sin(theta),
        col = "red", lty = 2, lwd = 1.5)
}

# One page per penalty
pdf(file.path(fig_dir, "dendrograms_application.pdf"), width = 14, height = 8)
for (k in seq_along(penalties)) plot_linear(trees[[k]], penalties[[k]])
invisible(dev.off())

pdf(file.path(fig_dir, "dendrograms_application_logscale.pdf"), width = 14, height = 8)
for (k in seq_along(penalties)) plot_logscale(trees[[k]], penalties[[k]])
invisible(dev.off())

pdf(file.path(fig_dir, "dendrograms_application_circular.pdf"), width = 12, height = 12)
for (k in seq_along(penalties)) plot_circular(trees[[k]], penalties[[k]])
invisible(dev.off())

# Figures 3, 4 and 5 of the paper (one file per penalty)
for (k in seq_along(penalties)) {
  pdf(file.path(fig_dir, sprintf("dendrograms_application_circular-%d.pdf", k)),
      width = 12, height = 12)
  plot_circular(trees[[k]], penalties[[k]])
  invisible(dev.off())
}
cat("Dendrograms saved to: results/application/figures/\n\n")

cat("=== Session information ===\n")
print(sessionInfo())

sink()
close(log_con)
