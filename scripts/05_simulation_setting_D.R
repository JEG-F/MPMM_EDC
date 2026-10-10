###############################################################################
#  05_simulation_setting_D.R
#
# Section 3, Setting D (spacing halved from one model to the next) -- Table 9
# of the paper.
# (Original script: simulation_study4.R.)
#
# A = {0,1}, p = 6, o = 5 (192 elements); 10 parts in every model:
#   Model 1 -> Equation (25): P(0|L_k) = 0.25   + (k-1) 0.05
#   Model 2 -> Equation (26): P(0|L_k) = 0.375  + (k-1) 0.025
#   Model 3 -> Equation (27): P(0|L_k) = 0.4375 + (k-1) 0.0125
# The partition of each model is drawn once (balanced random allocation with
# seeds 1001, 1002, 1003) and is the same in every replication.
#
# Usage (from the repository root):
#   Rscript scripts/05_simulation_setting_D.R               # as in the paper, 1 core
#   Rscript scripts/05_simulation_setting_D.R --cores=8     # same results, 8 cores
#   Rscript scripts/05_simulation_setting_D.R --cores=8 --models=1 --sizes=1000,9000
#   Rscript scripts/05_simulation_setting_D.R --reps=20     # quick test (other numbers)
#
# Output:  results/simulation/reruns/setting_D_errors_seed171266[...].csv, compared
#         cell by cell with results/simulation/setting_D_errors_seed171266.csv
#         (the totals used in the paper).  Table 9 itself is written by
#         scripts/06_simulation_tables.R.
#############################################################################

root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else normalizePath(".")
})
if (!file.exists(file.path(root, "R", "edc_partition.R")))
  stop("Run this script from the repository root: Rscript scripts/05_simulation_setting_D.R")

source(file.path(root, "R", "simulate_markov.R"))
source(file.path(root, "R", "edc_partition.R"))
source(file.path(root, "R", "simulation_tools.R"))

# --- Model definitions (as in simulation_study4.R) ---------------------------

m1_p0 <- 0.2500 + (seq_len(10) - 1) * 0.0500
m2_p0 <- 0.3750 + (seq_len(10) - 1) * 0.0250
m3_p0 <- 0.4375 + (seq_len(10) - 1) * 0.0125

models <- list(
  list(
    name        = "Model 1 (10 parts)",
    paper       = "Equation (25)",
    num_parts   = 10L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 10L,
                                     sim_partition_seeds[1]),
    trans_probs = make_trans_probs(m1_p0)
  ),
  list(
    name        = "Model 2 (10 parts)",
    paper       = "Equation (26)",
    num_parts   = 10L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 10L,
                                     sim_partition_seeds[2]),
    trans_probs = make_trans_probs(m2_p0)
  ),
  list(
    name        = "Model 3 (10 parts)",
    paper       = "Equation (27)",
    num_parts   = 10L,
    partition   = generate_partition(sim_p, sim_order_o, sim_alphabet, 10L,
                                     sim_partition_seeds[3]),
    trans_probs = make_trans_probs(m3_p0)
  )
)

run_simulation_study("D", "spacing halved from one model to the next",
                     models, root)
