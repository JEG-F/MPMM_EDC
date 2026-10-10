##################################################################################
#                        simulate_markov.R
#
# Simulate collections of independent discrete-time order-o Markov chains
# over a common finite alphabet, with a prescribed minimal partition
# structure (Definition 2-iii in the paper).
#
# The user specifies which (process, state) pairs share transition
# probabilities via a partition assignment, or lets the script generate
# a random partition.
#
# Usage:
#   Rscript simulate_markov.R <output_csv> <p> <N> <order_o> <A_size> <num_parts> [seed]
#
# Or from an R session:
#   source("simulate_markov.R")
#   sim <- simulate_markov_collection(p=4, N=5000, order_o=2,
#                                     alphabet=c("0","1"), num_parts=3)
#   write.csv(sim$data, "simulated.csv", row.names=FALSE)
##############################################################################

# --- Core simulation ---------------------------------------------------------

  simulate_markov_collection <- function(p, N, order_o, alphabet,
                                       partition = NULL,
                                       num_parts = NULL,
                                       trans_probs = NULL,
                                       min_prob = 0.05,
                                       min_sep = 0.15,
                                       max_attempts = 1000,
                                       seed = NULL) {
  # Simulate p independent order-o Markov chains on 'alphabet', each of
  # length N, with transition probabilities determined by a partition of
  # {1,...,p} x A^o.
  #
  # Args:
  #   p          : number of processes
  #   N          : sample size (length of each sequence)
  #   order_o    : Markov chain order
  #   alphabet   : character vector of symbols
  #   partition  : (optional) named list mapping part labels to vectors of
  #                element names "(j,s)" -- e.g. list(L1=c("(1,00)","(2,01)"), ...)
  #                If NULL, a random partition is generated.
  #   num_parts  :(optional) number of parts for the random partition.
  #                Ignored if 'partition' is given.
  #   trans_probs: (optional) matrix (num_parts x |A|) of transition
  #                probabilities, rows sum to 1. Row k is used for
  #                all (j,s) in the k-th part.
  #                If NULL, random transition probabilities are generated.
  #   min_prob   : minimum probability per symbol when generating random
  #                transition probabilities (ensures all symbols can occur).
  #   min_sep    : minimum total-variation distance between the transition
  #                probability vectors of any two distinct parts.
  #                TV(P,Q) = 0.5 * sum |P(a) - Q(a)|, ranges in [0, 1].
  #                Default 0.15 ensures parts are well-separated.
  #   max_attempts:  maximum number of sampling attempts to achieve min_sep.
  #   seed       :  random seed for reproducibility.
  #
  # Returns: a list with
  #   data       : data.frame (N rows x p columns) ready for write.csv
  #   partition  : the true partition (list of data.frames with sample, state)
  #   trans_probs: matrix of true transition probabilities per part
  #   part_assignment:  named vector mapping "(j,s)" -> part index

  if (!is.null(seed)) set.seed(seed)

  A_size <- length(alphabet)

  # All states in A^o
  state_grid <- expand.grid(rep(list(alphabet), order_o), stringsAsFactors = FALSE)
  states <- apply(state_grid, 1, paste0, collapse = "")
  num_states <- length(states)

  # All elements of {1,...,p} x A^o
  elements <- expand.grid(sample = seq_len(p), state = states,
                          stringsAsFactors = FALSE)
  elem_labels <- paste0("(", elements$sample, ",", elements$state, ")")
  num_elements <- nrow(elements)

  # -- Build or validate partition ---
  if (is.null(partition)) {
    if (is.null(num_parts)) {
      num_parts <- min(num_elements, max(2, ceiling(num_elements / 4)))
    }
    num_parts <- min(num_parts, num_elements)

    # Random partition:  assign each element to a random part,
    # ensuring every part gets at least one element.
    assignment <- sample(seq_len(num_parts), num_elements, replace = TRUE)
    # Fix: guarantee all parts are represented
    missing <- setdiff(seq_len(num_parts), unique(assignment))
    if (length(missing) > 0) {
      fix_idx <- sample(seq_len(num_elements), length(missing))
      assignment[fix_idx] <- missing
    }
  } else {
    # Parse provided partition
    num_parts <- length(partition)
    assignment <- integer(num_elements)
    for (k in seq_len(num_parts)) {
      for (el in partition[[k]]) {
        idx <- which(elem_labels == el)
        if (length(idx) == 0)
          stop("Element '", el, "' in partition not found among (j,s) pairs.")
        assignment[idx] <- k
      }
    }
    if (any(assignment == 0))
      stop("Some (j,s) elements are not assigned to any part.")
  }

  names(assignment) <- elem_labels

  # --- Build or validate transition probabilities ---
  if (is.null(trans_probs)) {
    # Generate random transition probabilities with minimum total-variation
    # distance min_sep between every pair of parts.
    # TV(P, Q) = 0.5 * sum_a |P(a) - Q(a)|
    trans_probs <- matrix(nrow = num_parts, ncol = A_size)
    colnames(trans_probs) <- alphabet

    draw_one <- function() {
      raw <- runif(A_size, min = min_prob, max = 1)
      raw / sum(raw)
    }

    tv_dist <- function(p_vec, q_vec) 0.5 * sum(abs(p_vec - q_vec))

    trans_probs[1, ] <- draw_one()
    for (k in seq(2, num_parts)) {
      accepted <- FALSE
      for (att in seq_len(max_attempts)) {
        candidate <- draw_one()
        seps <- vapply(seq_len(k - 1), function(i)
          tv_dist(candidate, trans_probs[i, ]), numeric(1))
        if (all(seps >= min_sep)) {
          accepted <- TRUE
          break
        }
      }
      if (!accepted)
        warning("Part ", k, ": could not achieve min_sep = ", min_sep,
                " after ", max_attempts, " attempts. Best separation: ",
                round(max(seps), 4),
                ". Consider reducing min_sep or num_parts.")
      trans_probs[k, ] <- candidate
    }
  } else {
    if (!is.matrix(trans_probs) || nrow(trans_probs) != num_parts ||
        ncol(trans_probs) != A_size)
      stop("trans_probs must be a matrix with ", num_parts, " rows and ",
           A_size, " columns.")
    if (any(abs(rowSums(trans_probs) - 1) > 1e-10))
      stop("Rows of trans_probs must sum to 1.")
  }

  # Build lookup: for each (j, s) -> transition probability vector
  # prob_lookup[[j]][[s]] = probability vector of length A_size
  prob_lookup <- vector("list", p)
  for (j in seq_len(p)) {
    prob_lookup[[j]] <- list()
    for (s in states) {
      lab <- paste0("(", j, ",", s, ")")
      k <- assignment[lab]
      prob_lookup[[j]][[s]] <- trans_probs[k, ]
    }
  }

  # --- Simulate each chain ---
  data_mat <- matrix("", nrow = N, ncol = p)
  colnames(data_mat) <- paste0("S", seq_len(p))

  for (j in seq_len(p)) {
    seq_j <- character(N)

    # Initial segment: draw first 'order_o' symbols uniformly
    seq_j[seq_len(order_o)] <- sample(alphabet, order_o, replace = TRUE)

    for (t in (order_o + 1):N) {
      s <- paste0(seq_j[(t - order_o):(t - 1)], collapse = "")
      probs <- prob_lookup[[j]][[s]]
      seq_j[t] <- sample(alphabet, 1, prob = probs)
    }

    data_mat[, j] <- seq_j
  }

  # --- Build partition as list of data.frames (same format as edc_partition) ---
  part_list <- vector("list", num_parts)
  for (k in seq_len(num_parts)) {
    idx <- which(assignment == k)
    part_list[[k]] <- elements[idx, , drop = FALSE]
    rownames(part_list[[k]]) <- NULL
  }

  list(data            = as.data.frame(data_mat, stringsAsFactors = FALSE),
       partition       = part_list,
       trans_probs     = trans_probs,
       part_assignment = assignment,
       alphabet        = alphabet,
       states          = states,
       order_o         = order_o)
}

# --- Display -----------------------------------------------------------------

print_simulation_summary <- function(sim) {
  cat("=== Simulation Summary ===\n")
  cat("Processes (p):", ncol(sim$data), "\n")
  cat("Sample size (N):", nrow(sim$data), "\n")
  cat("Order (o):", sim$order_o, "\n")
  cat("Alphabet:", paste(sim$alphabet, collapse = ", "), "\n")
  cat("States |A^o|:", length(sim$states), "\n")
  cat("True number of parts:", length(sim$partition), "\n\n")

  for (k in seq_along(sim$partition)) {
    members <- sim$partition[[k]]
    cat("--- Part", k, "---\n")
    cat("  P(a|L):", paste(paste0(sim$alphabet, "=",
                                   round(sim$trans_probs[k, ], 4)),
                            collapse = "  "), "\n")
    cat("  Elements:", nrow(members), " ->")
    show_n <- min(nrow(members), 8)
    for (i in seq_len(show_n)) {
      cat(" (", members$sample[i], ",", members$state[i], ")")
    }
    if (nrow(members) > show_n)
      cat(" ... [", nrow(members) - show_n, " more]")
    cat("\n")
  }
}

# --- CLI interface -----------------------------------------------------------

if (!interactive() && grepl("simulate_markov\\.R$", commandArgs()[4], ignore.case = TRUE)) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 6) {
    cat("Usage: Rscript simulate_markov.R <output_csv> <p> <N> <order_o>",
        "<A_size> <num_parts> [seed]\n\n")
    cat("  output_csv : path for the generated CSV\n")
    cat("  p          : number of processes (columns)\n")
    cat("  N          : sample size (rows)\n")
    cat("  order_o    : Markov chain order\n")
    cat("  A_size     : alphabet size (uses 0, 1, ..., A_size-1)\n")
    cat("  num_parts  : number of parts in the true partition\n")
    cat("  seed       : (optional) random seed\n")
    quit(status = 1)
  }

  output_csv <- args[1]
  p_val      <- as.integer(args[2])
  N_val      <- as.integer(args[3])
  order_o    <- as.integer(args[4])
  A_size     <- as.integer(args[5])
  num_parts  <- as.integer(args[6])
  seed_val   <- if (length(args) >= 7) as.integer(args[7]) else NULL

  alphabet <- as.character(seq(0, A_size - 1))

  sim <- simulate_markov_collection(p = p_val, N = N_val, order_o = order_o,
                                    alphabet = alphabet, num_parts = num_parts,
                                    seed = seed_val)

  write.csv(sim$data, output_csv, row.names = FALSE)
  cat("CSV written to:", output_csv, "\n\n")
  print_simulation_summary(sim)

  # Save the true partition to a file for reference
  truth_file <- sub("\\.csv$", "_truth.txt", output_csv)
  sink(truth_file)
  print_simulation_summary(sim)
  cat("\n=== Full partition assignment ===\n")
  for (k in seq_along(sim$partition)) {
    members <- sim$partition[[k]]
    cat("Part", k, ":\n")
    for (i in seq_len(nrow(members))) {
      cat("  (", members$sample[i], ",", members$state[i], ")\n")
    }
  }
  sink()
  cat("True partition saved to:", truth_file, "\n")
}
