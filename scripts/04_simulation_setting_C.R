###############################################################################
# 04_simulation_setting_C.R
#
#  Section 3, Setting C (well-separated parts plus a cluster of close part)
# ----- Table 8 of the paper.
# (Original script: simulation_study3.R.)
#
# A = {0,1}, p = 6, o = 5 (192 elements):
#   Model 1 -> not reported in the paper: 5 parts, P(0|L) = 0.15, 0.30, 0.50,
#              0.75, 0.80 (kept so that the stored result file is reproduced
#              in full; skip it with --models=2,3)
#   Model 2 -> Equations (18)-(20): 10 parts,
#              P(0|L_k) = 0.10 + (k-1) 0.09, k = 1..6, and a close cluster
#              P(0|L_{6+j}) = 0.775 + (j - 2.5) 0.0225, j = 1..4
#              (= 0.74125, 0.76375, 0.78625, 0.80875; printed with four
#              decimals in Equations (19)-(20))
#    Model 3 -> Equations (21)-(24): 15 parts,
#              P(0|L_k) = 0.08 + (k-1) 0.06, k = 1..9, and a close cluster
#              P(0|L_{9+j}) = 0.77 + (j - 3.5) 0.012, j = 1..6
#              (= 0.740, 0.752, 0.764, 0.776, 0.788, 0.800)
# The partitions are those of Setting B (same seeds and sizes; Table 7).
#
# Usage (from the repository root):
#   Rscript scripts/04_simulation_setting_C.R               # as in the paper, 1 core
#   Rscript scripts/04_simulation_setting_C.R --cores=8     # same results, 8 cores
#   Rscript scripts/04_simulation_setting_C.R --cores=8 --models=2,3
#   Rscript scripts/04_simulation_setting_C.R --reps=20     # quick test (other numbers)
#
#  Output: results/simulation/reruns/setting_C_errors_seed171266[...].csv, compared
#         cell by cell with results/simulation/setting_C_errors_seed171266.csv
#         (the totals used in the paper).  Table 8 itself is written by
#         scripts/06_simulation_tables.R.
###############################################################################

root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else normalizePath(".")
})
if (!file.exists(file.path(root, "R", "edc_partition.R")))
  stop("Run this script from the repository root: Rscript scripts/04_simulation_setting_C.R")

source(file.path(root, "R", "simulate_markov.R"))
source(file.path(root, "R", "edc_partition.R"))
source(file.path(root, "R", "simulation_tools.R"))

# --- Model definitions (as in simulation_study3.R) ---------------------------

# Close cluster for Model 2: 4 parts centered at 0.775 with spacing 0.0225
m2_close <- 0.775 + (seq_len(4) - 2.5) * 0.0225

# Close cluster for Model 3: 6 parts centered at 0.77 with spacing 0.012
m3_close <- 0.77 + (seq_len(6) - 3.5) * 0.012

models <- list(
  list(
    name        = "Model 1 (5 parts)",
    paper       = NULL,
    num_parts   = 5L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 5L,
                                     sim_partition_seeds[1]),
    trans_probs = make_trans_probs(c(0.15, 0.30, 0.50, 0.75, 0.80))
  ),
  list(
    name        = "Model 2 (10 parts)",
    paper       = "Equations (18)-(20)",
    num_parts   = 10L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 10L,
                                     sim_partition_seeds[2]),
    trans_probs = make_trans_probs(c(0.10, 0.19, 0.28, 0.37, 0.46, 0.55,
                                     m2_close))
  ),
  list(
    name        = "Model 3 (15 parts)",
    paper       = "Equations (21)-(24)",
    num_parts   = 15L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 15L,
                                     sim_partition_seeds[3]),
    trans_probs = make_trans_probs(c(0.08, 0.14, 0.20, 0.26, 0.32, 0.38,
                                     0.44, 0.50, 0.56, m3_close))
  )
)

run_simulation_study("C", "well-separated parts plus a cluster of close parts",
                     models, root)
