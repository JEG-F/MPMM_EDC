###############################################################################
# 00_download_data.R  (optional)
#
# Rebuilds data/dados_yahoo.RDATA from Yahoo Finance with the package
# 'tidyquant' and compares the new download with the stored copy.  The
# analyses always use the stored copy (downloaded in April 2026), because
# Yahoo Finance may revise or extend its historical series.  A download  made
# on 2026-10-08 gave exactly the same closing quotes (including the missing
# ones) and therefore the same binary series.
#
# Binary series (column 'sinal' of dados_yahoo).  For each currency, let
# y_1, ..., y_m be the closing quotes (US dollars per unit of the currency) of
# the m = 5182 trading days in the window, with a missing quote replaced by 0.
# Then
#     x_1 = 1   and   x_t = 1 if y_t >= y_{t-1},  x_t = 0 otherwise (t >= 2).
# When both quotes are available this is the rule of Section 4 of the paper.
# Around the 137 missing quotes (97 gaps, 0.44% of the currency-days) it gives x_t = 0
# on the first missing day of a gap and x_t = 1 on the remaining missing days
# and on the first quoted day after the gap.  See data/README.md.
#
# Usage (from the repository root): Rscript scripts/00_download_data.R
# Output: data/dados_yahoo_redownload.RDATA (the stored copy is not modified)
###############################################################################

root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else normalizePath(".")
})
if (!file.exists(file.path(root, "R", "edc_partition.R")))
  stop("Run this script from the repository root: Rscript scripts/00_download_data.R")
if (!requireNamespace("tidyquant", quietly = TRUE))
  stop("This script needs the package 'tidyquant': install.packages(\"tidyquant\")")

symbols   <- c("CNYUSD=X", "MXNUSD=X", "BRLUSD=X", "GBPUSD=X", "CADUSD=X", "EURUSD=X")
first_day <- as.Date("2006-03-24")
last_day  <- as.Date("2026-02-02")

raw <- tidyquant::tq_get(symbols, get = "stock.prices",
                         from = first_day, to = last_day + 1)
raw <- raw[raw$date >= first_day & raw$date <= last_day, ]
raw <- raw[order(raw$symbol, raw$date), ]

binary_signal <- function(close) {
  y <- close
  y[is.na(y)] <- 0
  c(1, as.numeric(diff(y) >= 0))
}
raw$sinal <- ave(raw$close, raw$symbol, FUN = binary_signal)

cat("Downloaded:", nrow(raw), "rows;",
    paste(names(table(raw$symbol)), table(raw$symbol), sep = ": ", collapse = ", "), "\n")

# ------ Comparison with the stored copy ----------------------------------------

e <- new.env()
load(file.path(root, "data", "dados_yahoo.RDATA"), envir = e)
old <- as.data.frame(e$dados_yahoo)
old <- old[order(old$symbol, old$date), ]
new <- as.data.frame(raw)

same_days <- identical(paste(old$symbol, old$date), paste(new$symbol, new$date))
cat("Same currencies and trading days as the stored copy:", same_days, "\n")
if (same_days) {
  cat("Closing quotes identical (missing ones included):",
      identical(old$close, new$close), "\n")
  cat("Binary series identical:",
      identical(as.numeric(old$sinal), as.numeric(new$sinal)), "\n")
} else {
  cat("The trading days differ; compare the two files before using the new one.\n")
}

dados_yahoo <- raw
out <- file.path(root, "data", "dados_yahoo_redownload.RDATA")
save(dados_yahoo, file = out)
cat("New download saved to:", out, "\n")
