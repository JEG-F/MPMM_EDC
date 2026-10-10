###############################################################################
# edc_partition.R
#
# EDC based estimation of minimal partitions for collections of independent 
# discrete-time Markov chains (Multivariate Partition Markov Models), via the 
# EDC-based dissimilarity d.
#
# Reference: D. S. Santos, J. E. Garcia & V. A. Gonzalez-Lopez (2026),
#   "Consistent Identification of Partition Markov Structures in Exchange
#   Rate Dynamics" (submitted to Applied Stochastic Models in Business and
#   Industry).
#
#  Correspondence with the paper (numbering of the 2026-10-08 manuscript):
#   compute_lml_part()  LML(L) of Equation (11); summed over the parts it gives
#                       the LML of Equation (6)
#   compute_K_inv()     K^{-1} of Definition 3
#   edc_partition()     Algorithm 1 (greedy partition recovery via d); it stops
#                       as soon as the smallest d is > 1 (Theorem 1)
#   $lml_value          LML, Equation (6), of the estimated partition
#   $edc_value          EDC, Equation (7), of the estimated partition
#
#  Usage:
#   Rscript R/edc_partition.R <csv_file> <order_o> [alpha]
#
#   Or from an R session:
#    source("R/edc_partition.R")
#    result <- run_edc("data.csv", order_o = 2)
###############################################################################

# --- Helpers -----------------------------------------------------------------

safe_xlogx <- function(x) {
  # x * log(x) with convention 0 * log(0) = 0, vectorized
  result <- numeric(length(x))
  pos <- x > 0
  result[pos] <- x[pos] * log(x[pos])
  result
}

compute_lml_part <- function(counts_a, count_total) {
  # LML contribution of a single part: sum_a N(L,a) * ln(N(L,a) / N(L))
  if (count_total == 0) return(0)
  sum(safe_xlogx(counts_a) - counts_a * log(count_total))
}

compute_K_inv <- function(A_size, num_parts, total_n, gamma_fn, w_fn, alpha) {
  # K^{-1} = (gamma((|A|-1)*|L|/alpha) - gamma((|A|-1)*(|L|-1)/alpha)) * w(sum n_j)
  (gamma_fn((A_size - 1) * num_parts / alpha) -
     gamma_fn((A_size - 1) * (num_parts - 1) / alpha)) * w_fn(total_n)
}

# --- Data processing ---------------------------------------------------------

compute_counts <- function(data, order_o) {
  # Extract N^j(s, a) counts from a data matrix (N rows x p columns).
  # Each column is a sample sequence of discrete symbols.
  data <- as.matrix(data)
  storage.mode(data) <- "character"
  p <- ncol(data)
  n <- nrow(data)
  alphabet <- sort(unique(as.vector(data)))
  A_size <- length(alphabet)

  # All A^o states
  state_grid <- expand.grid(rep(list(alphabet), order_o), stringsAsFactors = FALSE)
  states <- apply(state_grid, 1, paste0, collapse = "")

  # counts[j, s, a]
  counts <- array(0, dim = c(p, length(states), A_size),
                  dimnames = list(sample = seq_len(p),
                                  state  = states,
                                  symbol = alphabet))

  for (j in seq_len(p)) {
    seq_j <- data[, j]
    for (t in (order_o + 1):n) {
      s <- paste0(seq_j[(t - order_o):(t - 1)], collapse = "")
      a <- seq_j[t]
      counts[j, s, a] <- counts[j, s, a] + 1L
    }
  }

  state_counts <- apply(counts, c(1, 2), sum)   # N^j(s)

  list(alphabet     = alphabet,
       states       = states,
       p            = p,
       n            = rep(n, p),
       A_size       = A_size,
       counts       = counts,
       state_counts = state_counts)
}

# -- Partition initialisation ----------------------------------------------

initialize_partition <- function(counts_info) {
  # Finest partition: each (j, s) with N^j(s) > 0 is a singleton.
  p      <- counts_info$p
  states <- counts_info$states
  A_size <- counts_info$A_size

  parts           <- list()
  part_counts_a   <- matrix(nrow = 0, ncol = A_size)
  colnames(part_counts_a) <- counts_info$alphabet
  part_counts_total <- numeric(0)
  idx <- 0L

  for (j in seq_len(p)) {
    for (si in seq_along(states)) {
      s <- states[si]
      if (counts_info$state_counts[j, s] > 0) {
        idx <- idx + 1L
        parts[[idx]] <- data.frame(sample = j, state = s,
                                   stringsAsFactors = FALSE)
        part_counts_a   <- rbind(part_counts_a, counts_info$counts[j, s, ])
        part_counts_total <- c(part_counts_total,
                               counts_info$state_counts[j, s])
      }
    }
  }

  list(parts             = parts,
       part_counts_a     = part_counts_a,
       part_counts_total = part_counts_total)
}

# --Main algorithm -------------------------------------------------------

edc_partition <- function(counts_info,
                          gamma_fn = identity,
                          w_fn     = function(n) log(log(n)),
                          alpha    = 2) {

  init <- initialize_partition(counts_info)
  parts             <- init$parts
  part_counts_a     <- init$part_counts_a
  part_counts_total <- init$part_counts_total

  A_size  <- counts_info$A_size
  total_n <- sum(counts_info$n)

  if (total_n <= exp(1)) {
    warning("Total sample size is very small; w = log(log(n)) may be unstable.")
  }

  num_parts <- length(parts)
  merge_history <- list()

  #  Build initial delta matrix (K-independent raw divergences) ---
  # delta[v,m] = LML_v + LML_m - LML_{v union m}
  # Pre-compute per-part LML values
  lml_parts <- vapply(seq_len(num_parts), function(v)
    compute_lml_part(part_counts_a[v, ], part_counts_total[v]),
    numeric(1))

  delta <- matrix(Inf, nrow = num_parts, ncol = num_parts)
  if (num_parts > 1) {
    for (v in seq_len(num_parts - 1)) {
      for (m in (v + 1):num_parts) {
        ca_vm <- part_counts_a[v, ] + part_counts_a[m, ]
        ct_vm <- part_counts_total[v] + part_counts_total[m]
        lml_vm <- compute_lml_part(ca_vm, ct_vm)
        delta[v, m] <- lml_parts[v] + lml_parts[m] - lml_vm
        delta[m, v] <- delta[v, m]
      }
    }
  }

  # ---Greedy agglomerative loop ---
  step <- 0L
  repeat {
    num_parts <- length(parts)
    if (num_parts <= 1) break

    K_inv <- compute_K_inv(A_size, num_parts, total_n, gamma_fn, w_fn, alpha)

    # Find pair with minimum delta (=> minimum d since K_inv > 0 constant)
    min_idx <- which.min(delta)
    min_row <- ((min_idx - 1L) %% num_parts) + 1L
    min_col <- ((min_idx - 1L) %/% num_parts) + 1L
    min_d   <- delta[min_row, min_col] / K_inv

    if (min_d > 1) break

    step <- step + 1L
    v <- min(min_row, min_col)
    m <- max(min_row, min_col)

    merge_history[[step]] <- list(step = step, part_v = v, part_m = m,
                                  d_value = min_d)

    # Merge part m into part v
    parts[[v]] <- rbind(parts[[v]], parts[[m]])
    part_counts_a[v, ] <- part_counts_a[v, ] + part_counts_a[m, ]
    part_counts_total[v] <- part_counts_total[v] + part_counts_total[m]

    # Remove part m
    parts[[m]] <- NULL
    part_counts_a     <- part_counts_a[-m, , drop = FALSE]
    part_counts_total <- part_counts_total[-m]
    delta             <- delta[-m, -m, drop = FALSE]
    lml_parts         <- lml_parts[-m]

    # Update LML for merged part v
    lml_parts[v] <- compute_lml_part(part_counts_a[v, ], part_counts_total[v])

    # Recompute delta row/col for merged part v
    new_num <- length(parts)
    for (k in seq_len(new_num)) {
      if (k == v) {
        delta[v, k] <- Inf
        next
      }
      ca_vk <- part_counts_a[v, ] + part_counts_a[k, ]
      ct_vk <- part_counts_total[v] + part_counts_total[k]
      lml_vk <- compute_lml_part(ca_vk, ct_vk)
      delta[v, k] <- lml_parts[v] + lml_parts[k] - lml_vk
      delta[k, v] <- delta[v, k]
    }
  }

  # ---Final outputs------
  num_parts <- length(parts)
  if (num_parts > 1) {
    K_inv_final <- compute_K_inv(A_size, num_parts, total_n,
                                 gamma_fn, w_fn, alpha)
    d_matrix <- delta / K_inv_final
    diag(d_matrix) <- NA
  } else {
    d_matrix <- matrix(NA, 1, 1)
  }

  transition_probs <- part_counts_a / part_counts_total
  colnames(transition_probs) <- counts_info$alphabet

  lml_total <- sum(lml_parts)
  penalty   <- gamma_fn((A_size - 1) * num_parts / alpha) * w_fn(total_n)
  edc_value <- lml_total - penalty

  list(partition        = parts,
       transition_probs = transition_probs,
       d_matrix         = d_matrix,
       edc_value        = edc_value,
       lml_value        = lml_total,
       merge_history    = merge_history)
}

# --Entry point -------------------------------------------------------

run_edc <- function(csv_file, order_o,
                    gamma_fn = identity,
                    w_fn     = function(n) log(log(n)),
                    alpha    = 2) {

  data <- read.csv(csv_file, header = TRUE, colClasses = "character")

  cat("=== EDC Partition Algorithm ===\n")
  cat("File:          ", csv_file, "\n")
  cat("Samples (p):   ", ncol(data), "\n")
  cat("Sample size (N):", nrow(data), "\n")
  cat("Markov order (o):", order_o, "\n")

  counts_info <- compute_counts(data, order_o)

  cat("Alphabet:       ", paste(counts_info$alphabet, collapse = ", "), "\n")
  cat("States |A^o|:   ", length(counts_info$states), "\n")
  cat("Active (j,s):   ", sum(counts_info$state_counts > 0), "\n\n")

  result <- edc_partition(counts_info, gamma_fn, w_fn, alpha)

  cat("=== Results ===\n")
  cat("Parts in estimated partition:", length(result$partition), "\n")
  cat("LML value:                   ", round(result$lml_value, 4), "\n")
  cat("EDC value:                   ", round(result$edc_value, 4), "\n")
  cat("Merges performed:            ", length(result$merge_history), "\n\n")

  for (v in seq_along(result$partition)) {
    cat("--- Part", v, "---\n")
    members <- result$partition[[v]]
    cat("  Elements:", nrow(members), "\n")
    cat("  P(a|L)  :", paste(paste0(counts_info$alphabet, "=",
                                     round(result$transition_probs[v, ], 4)),
                              collapse = "  "), "\n")
    cat("  Members :")
    show_n <- min(nrow(members), 10)
    for (i in seq_len(show_n)) {
      cat(" (", members$sample[i], ",", members$state[i], ")")
    }
    if (nrow(members) > show_n)
      cat(" ... [", nrow(members) - show_n, " more]")
    cat("\n")
  }

  if (length(result$partition) > 1) {
    cat("\n=== Dissimilarity matrix d (at termination, all > 1) ===\n")
    print(round(result$d_matrix, 4))
  }

  invisible(result)
}

# --Terminal interface -------------------------------------------------------

if (!interactive() && grepl("edc_partition\\.R$", commandArgs()[4], ignore.case = TRUE)) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    cat("Usage: Rscript R/edc_partition.R <csv_file> <order_o> [alpha (default 2)]\n")
    quit(status = 1)
  }
  csv_file <- args[1]
  order_o  <- as.integer(args[2])
  alpha    <- if (length(args) >= 3) as.numeric(args[3]) else 2
  run_edc(csv_file, order_o, alpha = alpha)
}
