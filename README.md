# MPMM_EDC

R code and data that reproduce the results of

> D. S. Santos, J. E. García and V. A. González-López (2026).
> *Consistent Identification of Partition Markov Structures in Exchange Rate Dynamics.*
> Submitted to *Applied Stochastic Models in Business and Industry*.

A Multivariate Partition Markov Model describes $p$ independent Markov chains
of order $o$ on an alphabet $A$ by the minimal partition of
$\lbrace 1,\dots,p\rbrace\times A^o$, in which each part groups the (process, state)
pairs that share the same transition probabilities. The paper derives a
dissimilarity $d$ between parts from the Efficient Determination Criterion
(EDC). Merging two parts increases the EDC if and only if $d \le 1$
(Theorem 1). Algorithm 1 of the paper uses this rule to estimate the
partition by greedy agglomeration. This repository implements it in base R.

## Quick start

Requirements: R ≥ 3.6.0 (tested with R 4.5.2) and the package `ape`, used for
the circular dendrograms. The optional download script also needs `tidyquant`.

```sh
git clone https://github.com/JEG-F/MPMM_EDC.git
cd MPMM_EDC
Rscript -e 'install.packages("ape", repos = "https://cloud.r-project.org")'  # if needed
Rscript run_all.R
```

`run_all.R` takes a few seconds and does three things:

1. runs the exchange-rate application (Tables 10–16, Figures 3–5);
2. builds Tables 4, 6, 8 and 9 from the stored results of the simulation
   studies;
3. compares every number with the value printed in the paper
   (`tests/check_against_paper.R`) and ends with `All checks passed.`

Re-running the Monte Carlo studies takes hours; see
[Simulation studies](#simulation-studies-section-3).

## Where each result comes from

| In the paper | Script | Output |
|---|---|---|
| Table 4 (Setting A) | `scripts/02_simulation_setting_A.R`, then `scripts/06_simulation_tables.R` | `results/simulation/setting_A_errors_seed171266.csv`, `results/simulation/tables/table_setting_A.tex` |
| Table 6 (Setting B) | `scripts/03_simulation_setting_B.R`, then `scripts/06_simulation_tables.R` | `results/simulation/setting_B_errors_seed171266.csv`, `results/simulation/tables/table_setting_B.tex` |
| Table 8 (Setting C) | `scripts/04_simulation_setting_C.R`, then `scripts/06_simulation_tables.R` | `results/simulation/setting_C_errors_seed171266.csv`, `results/simulation/tables/table_setting_C.tex` |
| Table 9 (Setting D) | `scripts/05_simulation_setting_D.R`, then `scripts/06_simulation_tables.R` | `results/simulation/setting_D_errors_seed171266.csv`, `results/simulation/tables/table_setting_D.tex` |
| Tables 3, 5, 7 (models of Settings A–C) | "Model Summary" printed by `scripts/02`–`04` | console |
| Table 10 | `scripts/01_application_exchange_rates.R` | `results/application/summary_LML_EDC.csv`, `results/application/tables/summary_LML_EDC.tex` |
| Tables 11–12 ($w_1$ and $w_2$) | `scripts/01_application_exchange_rates.R` | `results/application/partitions.csv`, `results/application/tables/partition_w1.tex`, `partition_w2.tex` |
| Table 13 | `scripts/01_application_exchange_rates.R` | `results/application/tables/states_per_currency_w1.tex` |
| Tables 14–15 ($w_3$) | `scripts/01_application_exchange_rates.R` | `results/application/partitions.csv`, `results/application/tables/partition_w3.tex` |
| Table 16 | `scripts/01_application_exchange_rates.R` | `results/application/tables/partition_decomposition_w1_w3.tex` |
| Figures 3, 4, 5 | `scripts/01_application_exchange_rates.R` | `results/application/figures/dendrograms_application_circular-1.pdf`, `-2.pdf`, `-3.pdf` |

Tables 1–2, Figures 1–2 and Algorithm 1 are definitions or illustrations and
involve no computation. Table and equation numbers refer to the manuscript of
8 October 2026.

## Repository layout

```
R/
  edc_partition.R        Algorithm 1: counts, LML, K^{-1}, greedy merging with d, EDC
  edc_dendrogram.R       full agglomeration tree (dendrograms of Figures 3-5)
  simulate_markov.R      simulation of p independent order-o chains with a given partition
  simulation_tools.R     error counting and main loop shared by the four settings
scripts/
  00_download_data.R     optional: downloads the data again from Yahoo Finance
  01_application_exchange_rates.R              Section 4
  02_simulation_setting_A.R ... 05_simulation_setting_D.R   Section 3
  06_simulation_tables.R Tables 4, 6, 8 and 9
data/                    exchange-rate data, see data/README.md
results/                 outputs; the simulation totals used in the paper are stored here
tests/                   check_against_paper.R and the values printed in the paper
run_all.R                runs the scripts in order
```

## Simulation studies (Section 3)

All settings use $A=\lbrace 0,1\rbrace$, $p=6$, $o=5$ (192 elements), $\gamma$ = identity,
$\alpha=2$, the three penalties of Table 1, $N = 1000, 2000, \dots, 9000$ and
1000 replications. For every model, penalty and $N$, the setting scripts count
the Type 1 and Type 2 errors of Table 2 over the replications.
`scripts/06_simulation_tables.R` turns these counts into the proportions of
Equation (12).

The totals used in the paper are stored in
`results/simulation/setting_{A,B,C,D}_errors_seed171266.csv`. To compute them
again:

```sh
Rscript scripts/02_simulation_setting_A.R --cores=8    # Setting A on 8 cores
Rscript run_all.R --simulations --cores=8              # the four settings
```

Each run writes `results/simulation/reruns/setting_X_errors_seed171266.csv`.
Runs that cover only some models or sample sizes, or use fewer replications,
add the suffix `_reps<r>_partial`. With 1000 replications, every run reports
how many cells are identical to the stored totals. Options:

| Option | Default | Meaning |
|---|---|---|
| `--cores=k` | 1 | spread the replications over k cores (Linux and macOS; Windows uses 1 core) |
| `--models=2,3` | all | run only some models of the setting |
| `--sizes=1000,9000` | 1000,…,9000 | run only some sample sizes |
| `--reps=r` | 1000 | number of replications; r < 1000 gives a quick test with different numbers |
| `--seed=s` | 171266 | master seed |

A subset of cells reproduces the stored totals exactly. For example,
`Rscript scripts/05_simulation_setting_D.R --cores=8 --models=3 --sizes=1000`
takes about two minutes. The number of cores does not change the results.

**Running time.** A complete setting has 81 cells of 1000 replications each.
On an Apple M4 Pro with R 4.5.2, it takes roughly 6–7 hours on one core and
about 55 minutes with 12 cores. The four settings take about 3.5 hours with
12 cores.

**Reproducibility details.**

- `set.seed(171266)` followed by `sample.int(.Machine$integer.max, 1000)`
  gives one seed per replication. Replication $r$ of every model, penalty and
  $N$ uses the $r$-th seed, so the three penalties are applied to the same
  samples.
- The true partition of model $k$ is drawn once with seed $1000+k$
  (balanced allocation, `generate_partition()`). It is the same in all
  replications.
- Each replication sets its own seed, so serial and parallel runs give the
  same totals.
- Identical numbers require R ≥ 3.6.0, because the default sampling algorithm
  of R changed in version 3.6.0.
- Check made on 8 October 2026: every model of every setting was re-run with
  these scripts for $N = 1000$ and $N = 9000$, with 1000 replications each
  and 12 cores. All 72 cells were identical to the stored totals.

**Notes on the models.**

- In the result files, the models are named as in the original scripts.
  - Setting A: Models 1–3 are Equations (13)–(15).
  - Setting B: Models 2 and 3 are Equations (16) and (17).
  - Setting C: Models 2 and 3 are Equations (18)–(20) and (21)–(24).
  - Setting D: Models 1–3 are Equations (25)–(27).
- Settings B and C also contain a 5-part Model 1 that the paper does not
  report. It is kept so that the stored files are reproduced in full; use
  `--models=2,3` to skip it.
- In Setting C, the close parts of Model 2 have
  $P(0|L_7),\dots,P(0|L_{10}) = 0.74125, 0.76375, 0.78625, 0.80875$, printed
  with four decimals in Equations (19)–(20). Those of Model 3 are
  $0.740, 0.752, \dots, 0.800$.

## Exchange-rate application (Section 4)

`scripts/01_application_exchange_rates.R` reads `data/dados_yahoo.RDATA` and
builds the binary series of the six currencies, 5,182 symbols each
([data/README.md](data/README.md) gives the coding rule, including the
treatment of missing quotes). It then runs Algorithm 1 with $o=5$,
$\gamma$ = identity, $\alpha=2$, $n = \sum_j n_j = 31{,}092$ and the penalties
$w_1=\ln(n)$, $w_2=[\ln(\ln(n))]^2$ and $w_3=\ln(\ln(n))$. It prints the
estimated partitions and writes the tables and figures listed above
(Table 10):

| $w(\cdot)$ | LML | EDC | parts |
|---|---:|---:|---:|
| $\ln(\cdot)$ | −21,380.92 | −21,411.95 | 6 |
| $[\ln(\ln(\cdot))]^2$ | −21,380.92 | −21,397.29 | 6 |
| $\ln(\ln(\cdot))$ | −21,375.79 | −21,386.31 | 9 |

The dendrograms continue the merges of Algorithm 1 past $d=1$ until a single
part remains. The radius shows, on a $\ln(d+1)$ scale, the value of $d$ at
each merge, and the red dashed circle marks $d = 1$.

## Using the code with other data

```r
source("R/edc_partition.R")
x <- read.csv("my_series.csv", colClasses = "character")  # one column per process
counts <- compute_counts(x, order_o = 2)
fit <- edc_partition(counts, gamma_fn = identity,
                     w_fn = function(n) log(log(n)), alpha = 2)
fit$partition          # list of parts, each a data.frame of (sample, state) pairs
fit$transition_probs   # P(a | L) for every part
fit$lml_value          # LML, Equation (6)
fit$edc_value          # EDC, Equation (7)
```

`compute_counts()` takes series of equal length, one column per process. The
alphabet is the set of symbols found in the data.

## Origin of the code

The scripts are those used to produce the results of the paper. File paths
were adapted to this layout and a few features were added:

| File | Original | Changes |
|---|---|---|
| `R/edc_partition.R` | `edc_partition.R` | also returns `lml_value`; comments |
| `R/simulate_markov.R` | `simulate_markov.R` | none |
| `R/edc_dendrogram.R` | `edc_full_tree()` in `run_application.R` | none (moved) |
| `scripts/01_application_exchange_rates.R` | `run_application.R` | LML taken from `lml_value`; LaTeX tables, CSV files and single-page figures added |
| `scripts/02`–`05` | `simulation_study1.R`–`simulation_study4.R` | shared code moved to `R/simulation_tools.R`; optional multicore runs; comparison with the stored totals |
| `scripts/06_simulation_tables.R` | (new) | builds Tables 4, 6, 8 and 9 in the layout of the paper |

## Citation

```bibtex
@article{SantosGarciaGonzalezLopez2026,
  author  = {Santos, Diego S. and Garc{\'\i}a, Jes{\'u}s E. and Gonz{\'a}lez-L{\'o}pez, V. A.},
  title   = {Consistent Identification of Partition {M}arkov Structures in Exchange Rate Dynamics},
  journal = {Applied Stochastic Models in Business and Industry},
  year    = {2026},
  note    = {Submitted}
}
```

## Licence

The code is released under the GNU General Public License, version 3 or (at
your option) any later version; see [LICENSE](LICENSE). The exchange-rate data
come from Yahoo Finance; see [data/README.md](data/README.md) for their terms
of use.

Questions: please open an issue in this repository or contact the
corresponding author of the paper.
