#############################################################################
# 03_simulation_setting_B.R
#
# Section 3, Setting B (uniform separation) -- Table 6 of the paper.
# (Original script: simulation_study2.R.)
#
# A = {0,1}, p = 6, o = 5 (192 elements); P(0|L_k) evenly spaced in [0, 1]:
#   Model 1 -> not reported in the paper: 5 parts, P(0|L) = 0.15, 0.30, 0.50,
#              0.70, 0.85 (kept so that the stored result file is reproduced
#              in full; skip it with --models=2,3)
#   Model 2 -> Equation (16): 10 parts, P(0|L_k) = 0.10 + (k-1) 0.09
#   Model 3 -> Equation (17): 15 parts, P(0|L_k) = 0.08 + (k-1) 0.06
# The partition of each model is drawn once (balanced random allocation with
# seeds 1001, 1002, 1003) and is the same in every replication (Table 5).
#
# Usage (from the repository root):
#   Rscript scripts/03_simulation_setting_B.R               # as in the paper, 1 core
#   Rscript scripts/03_simulation_setting_B.R --cores=8      # same results, 8 cores
#   Rscript scripts/03_simulation_setting_B.R --cores=8 --models=2,3
#   Rscript scripts/03_simulation_setting_B.R --reps=20     # quick test (other numbers)
#
# Output: results/simulation/reruns/setting_B_errors_seed171266[...].csv, compared
#         cell by cell with results/simulation/setting_B_errors_seed171266.csv
#         (the totals used in the paper).  Table 6 itself is written by
#         scripts/06_simulation_tables.R.
###############################################################################

root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else normalizePath(".")
})
if (!file.exists(file.path(root, "R", "edc_partition.R")))
  stop("Run this script from the repository root: Rscript scripts/03_simulation_setting_B.R")

source(file.path(root, "R", "simulate_markov.R"))
source(file.path(root, "R", "edc_partition.R"))
source(file.path(root, "R", "simulation_tools.R"))

# --- Model definitions (as in simulation_study2.R) ---------------------------

models <- list(
  list(
    name        = "Model 1 (5 parts)",
    paper       = NULL,
    num_parts   = 5L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 5L,
                                     sim_partition_seeds[1]),
    trans_probs = make_trans_probs(c(0.15, 0.30, 0.50, 0.70, 0.85))
  ),
  list(
    name        = "Model 2 (10 parts)",
    paper       = "Equation (16)",
    num_parts   = 10L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 10L,
                                     sim_partition_seeds[2]),
    trans_probs = make_trans_probs(seq(0.10, 0.91, length.out = 10))
  ),
  list(
    name        = "Model 3 (15 parts)",
    paper       = "Equation (17)",
    num_parts   = 15L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 15L,
                                     sim_partition_seeds[3]),
    trans_probs = make_trans_probs(seq(0.08, 0.92, length.out = 15))
  )
)

run_simulation_study("B", "uniform separation", models, root)
