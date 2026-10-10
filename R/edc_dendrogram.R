###############################################################################
# edc_dendrogram.R
#
# Full agglomeration tree of Algorithm 1, used to draw the dendrograms of
# Figures 3-5 of the paper .
#
# edc_full_tree() repeats the greedy merges of edc_partition() (same Delta
# matrix, same tie-breaking by which.min, same K^{-1} ) but does not stop at
# d > 1: it keeps merging until a single part remains and records, for every
# merge, the value of d at which it took place.  Cutting the tree at d = 1
# therefore returns the partition estimated by edc_partition().
#
# The merge heights are made non-decreasing before building the "hclust"
# object (required for plotting).  In the exchange-rate application  only
# 3 of the 190 merges per penalty are affected, all below d = 1 and by less
# than 5e-5, so the plotted heights are the values of d.
#
# Requires: R/edc_partition.R (initialize_partition, compute_lml_part,
#           compute_K_inv).
############################################################################


edc_full_tree <- function(counts_info, col_names,
                          gamma_fn = identity,
                          w_fn     = function(n) log(log(n)),
                          alpha    = 2) {

  init <- initialize_partition(counts_info)
  part_counts_a     <- init$part_counts_a
  part_counts_total <- init$part_counts_total

  A_size  <- counts_info$A_size
  total_n <- sum(counts_info$n)
  n_init  <- length(init$parts)

  # Labels using currency names
  labels <- character(n_init)
  for (i in seq_len(n_init)) {
    m <- init$parts[[i]]
    labels[i] <- paste0("(", col_names[m$sample[1]], ",", m$state[1], ")")
  }

  if (n_init <= 1) {
    stop("Need at least 2 active (j,s) elements to build a dendrogram.")
  }

  lml_parts <- vapply(seq_len(n_init), function(v)
    compute_lml_part(part_counts_a[v, ], part_counts_total[v]),
    numeric(1))

  delta <- matrix(Inf, nrow = n_init, ncol = n_init)
  for (v in seq_len(n_init - 1)) {
    for (m2 in (v + 1):n_init) {
      ca_vm <- part_counts_a[v, ] + part_counts_a[m2, ]
      ct_vm <- part_counts_total[v] + part_counts_total[m2]
      lml_vm <- compute_lml_part(ca_vm, ct_vm)
      delta[v, m2] <- lml_parts[v] + lml_parts[m2] - lml_vm
      delta[m2, v] <- delta[v, m2]
    }
  }

  cluster_id <- -(seq_len(n_init))
  merge_mat <- matrix(0L, nrow = n_init - 1, ncol = 2)
  heights   <- numeric(n_init - 1)
  num_parts <- n_init

  for (step in seq_len(n_init - 1)) {
    K_inv <- compute_K_inv(A_size, num_parts, total_n, gamma_fn, w_fn, alpha)

    min_idx <- which.min(delta)
    min_row <- ((min_idx - 1L) %% num_parts) + 1L
    min_col <- ((min_idx - 1L) %/% num_parts) + 1L
    d_val   <- delta[min_row, min_col] / K_inv

    v <- min(min_row, min_col)
    m <- max(min_row, min_col)

    merge_mat[step, ] <- sort(c(cluster_id[v], cluster_id[m]))
    heights[step] <- d_val

    cluster_id[v] <- step
    cluster_id <- cluster_id[-m]

    part_counts_a[v, ] <- part_counts_a[v, ] + part_counts_a[m, ]
    part_counts_total[v] <- part_counts_total[v] + part_counts_total[m]
    part_counts_a     <- part_counts_a[-m, , drop = FALSE]
    part_counts_total <- part_counts_total[-m]
    delta             <- delta[-m, -m, drop = FALSE]
    lml_parts         <- lml_parts[-m]

    lml_parts[v] <- compute_lml_part(part_counts_a[v, ], part_counts_total[v])

    num_parts <- num_parts - 1L
    if (num_parts <= 1) break

    for (k in seq_len(num_parts)) {
      if (k == v) { delta[v, k] <- Inf; next }
      ca_vk <- part_counts_a[v, ] + part_counts_a[k, ]
      ct_vk <- part_counts_total[v] + part_counts_total[k]
      lml_vk <- compute_lml_part(ca_vk, ct_vk)
      delta[v, k] <- lml_parts[v] + lml_parts[k] - lml_vk
      delta[k, v] <- delta[v, k]
    }
  }

  for (i in seq_along(heights)[-1]) {
    if (heights[i] < heights[i - 1]) heights[i] <- heights[i - 1]
  }

  hc <- list(merge  = merge_mat,
             height = heights,
             order  = seq_len(n_init),
             labels = labels,
             method = "EDC")
  class(hc) <- "hclust"

  dend <- as.dendrogram(hc)
  hc$order <- order.dendrogram(dend)

  hc
}


