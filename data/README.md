# Data

Daily exchange rates of the six currencies against the US dollar, 2006-03-24 to
2026-02-02 (5,182 trading days per currency), used in Section 4 of the paper.

| File | Content |
|---|---|
| `dados_yahoo.RDATA` | Original data used in the paper (R object `dados_yahoo`) |
| `exchange_rates_2006-03-24_2026-02-02.csv` | The same data as CSV, for use outside R |
| `binary_series.csv` | Binary series analysed in the paper (input of the algorithm) |

## Source

Yahoo Finance (<https://finance.yahoo.com>), tickers `BRLUSD=X`, `CADUSD=X`,
`CNYUSD=X`, `EURUSD=X`, `GBPUSD=X` and `MXNUSD=X` (US dollars per unit of the
currency), downloaded with `tidyquant::tq_get()` in April 2026.

`scripts/00_download_data.R` downloads the series again and compares them with
the stored copy. A new download on 2026-10-08 gave exactly the same closing
quotes (missing ones included) and the same binary series. Yahoo Finance may
revise its historical data, so the analyses always use the stored copy.

## `dados_yahoo.RDATA`

`load("data/dados_yahoo.RDATA")` creates `dados_yahoo`, a tibble with
31,092 rows (6 currencies x 5,182 days) and 9 columns:

| Column | Description |
|---|---|
| `symbol` | Yahoo Finance ticker |
| `date` | trading day (`Date`) |
| `open`, `high`, `low`, `close` | daily quotes, US dollars per unit of the currency |
| `volume` | always 0 for exchange rates |
| `adjusted` | equal to `close` |
| `sinal` | binary series $x_t$ of the paper (1 = increase or maintenance, 0 = decrease) |

MD5 checksum: `d73027c65a4b5e9391f50b1f8900286f`.

  The CSV file holds the same table. It was written with
`write.csv(as.data.frame(dados_yahoo), "data/exchange_rates_2006-03-24_2026-02-02.csv", row.names = FALSE)`.

## Binary series

For each currency, let $y_1,\dots,y_m$ be the closing quotes of the
$m = 5182$ trading days, with a missing quote replaced by 0. Then

```math
x_1 = 1, \qquad x_t = \begin{cases} 1 & \text{if } y_t \ge y_{t-1},\\ 0 & \text{if } y_t < y_{t-1},\end{cases} \qquad t = 2,\dots,m .
```

When both quotes are available, this is the rule of Section 4 of the paper.
This rule reproduces the column `sinal` exactly for all 31,092 rows.

Some closing quotes are missing in Yahoo Finance: 137 of the 31,092
currency-days (0.44%), in 97 gaps. Around a gap the rule gives $x_t = 0$ on the first missing day,
and $x_t = 1$ on the other missing days and on the first quoted day after the
gap:

| Currency | Missing closing quotes | Gaps | $x_t = 1$ | $x_t = 0$ |
|---|---:|---:|---:|---:|
| BRL | 42 | 24 | 2,625 | 2,557 |
| CAD | 17 | 16 | 2,584 | 2,598 |
| CNY | 15 | 14 | 2,828 | 2,354 |
| EUR | 32 | 17 | 2,611 | 2,571 |
| GBP | 20 | 16 | 2,600 | 2,582 |
| MXN | 11 | 10 | 2,697 | 2,485 |

`binary_series.csv` is the 5,182 x 6 matrix of $x_t$, one column per currency
in the order BRL, EUR, CAD, GBP, MXN, CNY, which is also the process order
$j = 1,\dots,6$ used by the code. The estimated partitions do not depend on
this order. The file is written by `scripts/01_application_exchange_rates.R`.
The algorithm uses all $n_j = 5182$ symbols of each series, so
$n = \sum_j n_j = 31{,}092$ in the penalty $w(n)$.

## Terms of use

The quotes come from Yahoo Finance. They are included only so that the
results of the paper can be reproduced, and Yahoo's terms of use apply to
them. The GPL-3.0 licence of this repository covers the code, not these data.
