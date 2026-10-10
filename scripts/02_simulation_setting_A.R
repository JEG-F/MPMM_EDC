###############################################################################
# 02_simulation_setting_A.R
#
# Section 3, Setting A (weak separation) -- Table 4 of the paper.
# (Original script: simulation_study1.R.)
#
# A = {0,1}, p = 6, o = 5 (192 elements); every P(0|L_k) lies in [0.40, 0.60]:
#   Model 1 -> Equation (13):  5 parts, P(0|L_k) = 0.40 + (k-1) 0.05
#   Model 2 -> Equation (14): 10 parts, P(0|L_k) = 0.40 + (k-1)/9  0.20
#   Model 3 -> Equation (15): 15 parts, P(0|L_k) = 0.40 + (k-1)/14 0.20
# The partition of each model is drawn once (balanced random allocation with
# seeds 1001, 1002, 1003) and is the same in every replication (Table 3).
#
#  Usage (from the repository root):
#   Rscript scripts/02_simulation_setting_A.R               # as in the paper, 1 core
#   Rscript scripts/02_simulation_setting_A.R --cores=8     # same results, 8 cores
#   Rscript scripts/02_simulation_setting_A.R --cores=8 --models=1 --sizes=1000,9000
#                                                           # a subset of the cells
#   Rscript scripts/02_simulation_setting_A.R --reps=20     # quick test (other numbers)
#
#   Output: results/simulation/reruns/setting_A_errors_seed171266[...].csv, compared
#         cell by cell with results/simulation/setting_A_errors_seed171266.csv
#         (the totals used in the paper).  Table 4 itself is written by
#         scripts/06_simulation_tables.R.
###############################################################################

root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else normalizePath(".")
})
if (!file.exists(file.path(root, "R", "edc_partition.R")))
  stop("Run this script from the repository root: Rscript scripts/02_simulation_setting_A.R")

source(file.path(root, "R", "simulate_markov.R"))
source(file.path(root, "R", "edc_partition.R"))
source(file.path(root, "R", "simulation_tools.R"))

# ------ Model definitions (as in simulation_study1.R) ----------------------

models <- list(
  list(
    name        = "Model 1 (5 parts)",
    paper       = "Equation (13)",
    num_parts   = 5L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 5L,
                                     sim_partition_seeds[1]),
    trans_probs = make_trans_probs(seq(0.40, 0.60, length.out = 5))
  ),
  list(
    name        = "Model 2 (10 parts)",
    paper       = "Equation (14)",
    num_parts   = 10L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 10L,
                                     sim_partition_seeds[2]),
    trans_probs = make_trans_probs(seq(0.40, 0.60, length.out = 10))
  ),
  list(
    name        = "Model 3 (15 parts)",
    paper       = "Equation (15)",
    num_parts   = 15L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 15L,
                                     sim_partition_seeds[3]),
    trans_probs = make_trans_probs(seq(0.40, 0.60, length.out = 15))
  )
)

run_simulation_study("A", "weak separation, all P(0|L) in [0.40, 0.60]",
                     models, root)
