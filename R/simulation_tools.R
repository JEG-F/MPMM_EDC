#############################################################################
#     simulation_tools.R
#
# Shared code of the four Monte Carlo studies of Section 3 of the paper
# (Settings A-D, scripts/02-05).  Every setting uses A = {0,1}, p = 6, o = 5
# (the 192 elements of {1,...,6} x {0,1}^5), the three penalties of Table 1,
# gamma = identity, alpha = 2, N = 1000, 2000, ..., 9000 and 1000 replications.
#
# count_errors(), generate_partition() and make_trans_probs() are copied
# verbatim from the original study scripts (simulation_study1.R ...
# simulation_study4.R).  run_simulation_study() keeps the original main loop
# (models x penalties x N x replications, replication r seeded with
# all_seeds[r]).  The only addition is that the replications of each
# (model, penalty, N) cell can be spread over several cores (--cores=k):
# every replication re-seeds the random number generator inside
# simulate_markov_collection(), so serial and parallel runs give identical
# totals.
#
# Requires: R/simulate_markov.R and R/edc_partition.R.
###############################################################################

# --- Error counting ----------------------------------------------------------
# Type 1 = pair in different true parts but in the same estimated part
#         (incorrectly together); Type 2 = pair in the same true part but in
#         different estimated parts ( incorrectly separated). Table 2.

count_errors <- function(true_partition, recovered_partition) {
  true_lab <- character(0)
  for (k in seq_along(true_partition)) {
    members <- true_partition[[k]]
    for (i in seq_len(nrow(members))) {
      key <- paste0("(", members$sample[i], ",", members$state[i], ")")
      true_lab[key] <- k
    }
  }
  rec_lab <- character(0)
  for (k in seq_along(recovered_partition)) {
    members <- recovered_partition[[k]]
    for (i in seq_len(nrow(members))) {
      key <- paste0("(", members$sample[i], ",", members$state[i], ")")
      rec_lab[key] <- k
    }
  }

  common <- intersect(names(true_lab), names(rec_lab))
  n <- length(common)
  t_vec <- true_lab[common]
  r_vec <- rec_lab[common]

  type1 <- 0L
  type2 <- 0L

  if (n > 1) {
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        same_true <- (t_vec[i] == t_vec[j])
        same_rec  <- (r_vec[i] == r_vec[j])
        if (!same_true && same_rec) type1 <- type1 + 1L
        if (same_true && !same_rec)  type2 <- type2 + 1L
      }
    }
  }

  c(type1 = type1, type2 = type2)
}

#  Generate a fixed random partition of {1,...,p} x A^o ---------------
# Balanced allocation (part sizes differ by at most one), randomly shuffled
# with a fixed seed; the partition is the same in every replication.

generate_partition <- function(p, order_o, alphabet, num_parts, seed) {
  set.seed(seed)
  states <- apply(expand.grid(rep(list(alphabet), order_o),
                              stringsAsFactors = FALSE),
                  1, paste0, collapse = "")
  elements <- expand.grid(sample = seq_len(p), state = states,
                          stringsAsFactors = FALSE)
  n_elem <- nrow(elements)
  assignments <- rep(seq_len(num_parts), length.out = n_elem)
  assignments <- sample(assignments)
  partition <- list()
  for (k in seq_len(num_parts)) {
    idx <- which(assignments == k)
    members <- character(length(idx))
    for (i in seq_along(idx)) {
      members[i] <- paste0("(", elements$sample[idx[i]], ",",
                           elements$state[idx[i]], ")")
    }
    partition[[paste0("L", k)]] <- members
  }
  partition
}

# ----Transition probabilities for a binary alphabet -----------------------
# Row k of the result is (P(0|L_k), P(1|L_k)).

make_trans_probs <- function(p0_vec) {
  matrix(c(p0_vec, 1 - p0_vec), ncol = 2,
         dimnames = list(NULL, c("0", "1")))
}

# --- Common design of the four settings --------------------------------------

sim_p               <- 6L
sim_order_o         <- 5L
sim_alphabet        <- c("0", "1")
sim_partition_seeds <- c(1001L, 1002L, 1003L)   # one per model, as in the originals
sim_sample_sizes    <- c(1000L, 2000L, 3000L, 4000L, 5000L, 6000L, 7000L, 8000L, 9000L)
sim_n_rep           <- 1000L
sim_seed            <- 171266L                  # master seed of the paper

# Penalties w(n) of Table 1 (labels as in the stored result files)
sim_penalties <- list(
  list(short = "ln(n)",
       w_fn  = function(n) log(n)),
  list(short = "(ln(ln(n)))^2",
       w_fn  = function(n) log(log(n))^2),
  list(short = "ln(ln(n))",
       w_fn  = function(n) log(log(n)))
)

# --- Command-line options ----------------------------------------------------
#   --seed=171266   master seed (a bare integer is also accepted, as in the
#                   original scripts)
#   --cores=1       number of cores for the replications (fork-based, so
#                   Linux/macOS only; Windows runs serially)
#   --reps=1000     number of replications; a smaller value uses the first
#                    replication seeds of the full run (quick test)
#   --models=1,2,3  subset of models of the setting
#    --sizes=1000,9000  subset of sample sizes N

parse_sim_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- list(seed = sim_seed, cores = 1L, reps = sim_n_rep,
              models = NULL, sizes = NULL)
  for (a in args) {
    if (grepl("^[0-9]+$", a)) {
      opt$seed <- as.integer(a)
      next
    }
    kv <- regmatches(a, regexec("^--([a-z]+)=([0-9,]+)$", a))[[1]]
    if (length(kv) != 3) stop("Unrecognised argument: ", a)
    val <- as.integer(strsplit(kv[3], ",")[[1]])
    switch(kv[2],
           seed   = opt$seed   <- val[1],
           cores  = opt$cores  <- val[1],
           reps   = opt$reps   <- val[1],
           models = opt$models <- val,
           sizes  = opt$sizes  <- val,
           stop("Unknown option: --", kv[2]))
  }
  opt
}

#  One cell: all replications for a (model, penalty, N) ---------------
#  Same steps as the body of the replication loop of the original scripts;
#  returns the total number of Type 1 and Type 2 errors over the replications.

run_replications <- function(mod, pen, N, seeds, cores = 1L) {
  one_rep <- function(seed_r) {
    sim <- simulate_markov_collection(
      p        = sim_p,
      N        = N,
      order_o  = sim_order_o,
      alphabet = sim_alphabet,
      partition   = mod$partition,
      trans_probs = mod$trans_probs,
      seed     = seed_r
    )

    counts_info <- compute_counts(sim$data, sim_order_o)
    result <- edc_partition(counts_info, gamma_fn = identity,
                            w_fn = pen$w_fn, alpha = 2)

    count_errors(sim$partition, result$partition)
  }

  if (cores > 1L) {
    errs <- parallel::mclapply(seeds, one_rep, mc.cores = cores)
    failed <- vapply(errs, inherits, logical(1), what = "try-error")
    if (any(failed)) stop("A replication failed: ", errs[[which(failed)[1]]])
    errs <- do.call(rbind, errs)
  } else {
    errs <- t(vapply(seeds, one_rep, c(type1 = 0L, type2 = 0L)))
  }
  c(type1 = sum(errs[, "type1"]), type2 = sum(errs[, "type2"]))
}

#  Printing helpers ------------------------------------------------------

print_model_summary <- function(models) {
  cat("=== Model Summary ===\n\n")
  for (mod in models) {
    sizes <- vapply(mod$partition, length, integer(1))
    p0    <- mod$trans_probs[, 1]
    gaps  <- diff(sort(p0))
    cat(mod$name, "--", if (is.null(mod$paper)) "not reported in the paper"
                        else paste("paper:", mod$paper), "\n")
    cat("  Elements:    ", sum(sizes), "\n")
    cat("  Part sizes:  ", paste(sizes, collapse = ", "), "\n")
    cat("  P(0|L):      ", paste(round(p0, 4), collapse = ", "), "\n")
    cat("  Spacings:    ", paste(round(gaps, 4), collapse = ", "), "\n")
    cat("  Min spacing: ", round(min(gaps), 4), "\n\n")
  }
}

compare_with_reference <- function(results, ref_file) {
  if (!file.exists(ref_file)) {
    cat("No stored results to compare with:", ref_file, "\n")
    return(invisible(NA))
  }
  ref  <- read.csv(ref_file, stringsAsFactors = FALSE)
  key  <- function(d) paste(d$Model, d$Penalty, d$N, sep = " | ")
  m    <- match(key(results), key(ref))
  same <- !is.na(m) & results$Type1 == ref$Type1[m] & results$Type2 == ref$Type2[m]
  same[is.na(same)] <- FALSE
  cat(sprintf("\nComparison with the stored results (%s): %d of %d cells identical.\n",
              basename(ref_file), sum(same), length(same)))
  if (!all(same)) {
    cat("Cells that differ:\n")
    print(cbind(results[!same, ], Stored_Type1 = ref$Type1[m[!same]],
                Stored_Type2 = ref$Type2[m[!same]]), row.names = FALSE)
  }
  invisible(all(same))
}

# - Main ----------------------------------------------------------

run_simulation_study <- function(setting, title, models, root,
                                 opt = parse_sim_args()) {
  penalties    <- sim_penalties
  n_rep        <- opt$reps
  seed         <- opt$seed
  sel          <- if (is.null(opt$models)) seq_along(models) else opt$models
  sample_sizes <- if (is.null(opt$sizes)) sim_sample_sizes else as.integer(opt$sizes)
  if (!all(sel %in% seq_along(models))) stop("--models must be in 1..", length(models))
  full_run <- identical(sel, seq_along(models)) &&
    identical(sample_sizes, sim_sample_sizes) && n_rep == sim_n_rep
  cores <- max(1L, opt$cores)
  if (cores > 1L && .Platform$OS.type == "windows") {
    message("Forked parallelism is not available on Windows; running on 1 core.")
    cores <- 1L
  }

  print_model_summary(models)

  cat(sprintf("=== Simulation study: Setting %s (%s) ===\n", setting, title))
  cat("Seed:          ", seed, "\n")
  cat("Replications:  ", n_rep, "\n")
  cat("Sample sizes:  ", paste(sample_sizes, collapse = ", "), "\n")
  cat("Models:        ", paste(vapply(models[sel], `[[`, "", "name"),
                               collapse = "; "), "\n")
  cat("Penalties:     ", length(penalties), "\n")
  cat("Cores:         ", cores, "\n\n")

  # Replication seeds exactly as in the original scripts: 1000 seeds drawn
  # after set.seed(seed); a run with --reps < 1000 uses the first ones.
  set.seed(seed)
  all_seeds <- sample.int(.Machine$integer.max, max(n_rep, sim_n_rep))[seq_len(n_rep)]

  results <- data.frame(
    Model   = character(0),
    Penalty = character(0),
    N       = integer(0),
    Type1   = integer(0),
    Type2   = integer(0),
    stringsAsFactors = FALSE
  )

  total_combos <- length(sel) * length(penalties) * length(sample_sizes)
  combo <- 0L

  for (mod in models[sel]) {
    for (pen in penalties) {
      for (N in sample_sizes) {
        combo <- combo + 1L
        cat(sprintf("[%d/%d] %s | %-15s | N=%d ... ",
                    combo, total_combos, mod$name, pen$short, N))
        t0 <- Sys.time()

        errs <- run_replications(mod, pen, N, all_seeds, cores)

        cat(sprintf("Type1=%d  Type2=%d  (%.0f s)\n", errs[["type1"]], errs[["type2"]],
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))

        results <- rbind(results, data.frame(
          Model   = mod$name,
          Penalty = pen$short,
          N       = N,
          Type1   = errs[["type1"]],
          Type2   = errs[["type2"]],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  #  Output tables (as in the original scripts) -----------------------

  cat("\n\n")
  cat(strrep("=", 78), "\n")
  cat("  RESULTS: Total errors across", n_rep, "replications\n")
  cat(strrep("=", 78), "\n")

  for (mod in models[sel]) {
    cat("\n")
    cat(strrep("-", 78), "\n")
    cat("  ", mod$name, "\n")
    cat(strrep("-", 78), "\n")
    cat(sprintf("  %-15s %6s %10s %10s %10s %10s\n",
                "Penalty", "N", "Type1", "Type2", "Total", "Avg/rep"))
    cat(strrep("-", 78), "\n")

    sub <- results[results$Model == mod$name, ]
    for (i in seq_len(nrow(sub))) {
      row <- sub[i, ]
      total_err <- row$Type1 + row$Type2
      avg_err   <- total_err / n_rep
      cat(sprintf("  %-15s %6d %10d %10d %10d %10.3f\n",
                  row$Penalty, row$N, row$Type1, row$Type2,
                  total_err, avg_err))
    }
    cat(strrep("-", 78), "\n")
  }

  # Save and compare with the totals used in the paper --------------------

  sim_dir <- file.path(root, "results", "simulation")
  out_dir <- file.path(sim_dir, "reruns")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, sprintf("setting_%s_errors_seed%d%s.csv", setting, seed,
                                         if (full_run) "" else sprintf("_reps%d_partial", n_rep)))
  write.csv(results, out_file, row.names = FALSE)
  cat("\nResults saved to:", out_file, "\n")

  if (n_rep == sim_n_rep) {
    compare_with_reference(results, file.path(sim_dir, sprintf("setting_%s_errors_seed%d.csv",
                                                               setting, seed)))
  } else {
    cat("(Comparison with the stored results requires --reps=", sim_n_rep, ")\n", sep = "")
  }
  invisible(results)
}
