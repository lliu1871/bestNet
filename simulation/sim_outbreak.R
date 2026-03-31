library(ape)
library(TransPhylo)
library(igraph)
library(networkDynamic)
library(ndtv)
library(network)
library(rgexf)
library(phybase)

simulate_outbreak <- function(initial_infection_time = 0, infection_rate = 2.0, removal_rate = 2.0, latent_mean = 0.2, target_size = 50) {
  repeat {
    outbreak <- data.frame(
      id = integer(),
      parent_id = integer(),
      infection_time = numeric(),
      latent_period = numeric(),
      onsite_time = numeric(),
      infectious_period = numeric(),
      removal_time = numeric(),
      num_infections = integer(),
      stringsAsFactors = FALSE
    )
    # queue of individuals to simulate
    queue <- data.frame(
      id = 1,
      parent_id = NA,
      infection_time = initial_infection_time
    )

    next_id <- 2

    while (nrow(queue) > 0) {
      # take the first individual
      person <- queue[1, ]
      queue <- queue[-1, ]

      # Step 1: latent period
      latent_period <- rchisq(1, df = latent_mean)

      # Step 2: on-site time
      onsite_time <- person$infection_time + latent_period

      # Step 3: infectious period
      infectious_period <- rexp(1, rate = removal_rate)

      # Step 4: removal time
      removal_time <- onsite_time + infectious_period

      # Step 5: number of infections
      num_infections <- rpois(1, lambda = infection_rate * infectious_period)

      # Step 6: secondary infection times only if less than 2*target_size
      if (num_infections > 0 && next_id < 2 * target_size) {
        secondary_times <- onsite_time + runif(num_infections, min = 0, max = infectious_period)
        new_infections <- data.frame(
          id = next_id:(next_id + num_infections - 1),
          parent_id = person$id,
          infection_time = secondary_times
        )
        queue <- rbind(queue, new_infections)
        next_id <- next_id + num_infections
      }

      # save this individual
      outbreak <- rbind(
        outbreak,
        data.frame(
          id = person$id,
          parent_id = person$parent_id,
          infection_time = person$infection_time,
          latent_period = latent_period,
          onsite_time = onsite_time,
          infectious_period = infectious_period,
          removal_time = removal_time,
          num_infections = num_infections
        )
      )
    }

    # check outbreak size
    if (nrow(outbreak) >= target_size) {
      # sort by onsite_time (ascending) before returning
      outbreak <- outbreak[order(outbreak$onsite_time), ]
      # keep a copy of old ids so we can remap parent_id after reindexing
      old_id <- outbreak$id
      # reset row names to be sequential
      rownames(outbreak) <- NULL
      # replace id with row numbers (1..n) and remap parent_id accordingly
      outbreak$id <- seq_len(nrow(outbreak))
      parent_old <- outbreak$parent_id
      # create lookup from old id -> new id (names are characters)
      lookup <- setNames(outbreak$id, as.character(old_id))
      # map parent_old to new ids; keep NA where parent_old is NA or not found
      outbreak$parent_id <- ifelse(is.na(parent_old), NA_integer_, as.integer(lookup[as.character(parent_old)]))
      # warn if any non-NA parent_old could not be remapped
      missing_parents <- !is.na(parent_old) & is.na(outbreak$parent_id)
      if (any(missing_parents)) {
        warning(sprintf("%d parent_id(s) were not found in outbreak after reindexing and were set to NA.", sum(missing_parents)))
      }

      # reset infection_time
      outbreak$infection_time <- outbreak$infection_time - outbreak$onsite_time[1]
      outbreak$removal_time <- outbreak$removal_time - outbreak$onsite_time[1]
      outbreak$onsite_time <- outbreak$onsite_time - outbreak$onsite_time[1]
      return(outbreak[1:target_size, ])
    }
    # else, repeat simulation until size ≥ target_size
  }
}

# Find all infectors (ancestors) of a given id in the outbreak
# outbreak: data.frame with columns 'id' and 'parent_id'
# id: integer id
# include_self: if TRUE, include the id itself as the first element
# returns: integer vector of ancestor ids ordered from immediate parent up to the root
find_infectors <- function(outbreak, id, include_self = TRUE) {
  if (!is.data.frame(outbreak) || !all(c("id", "parent_id") %in% names(outbreak))) stop("outbreak must be a data.frame with 'id' and 'parent_id'")
  id <- as.integer(id)
  if (!(id %in% outbreak$id)) {
    warning("id not found in outbreak")
    return(integer(0))
  }
  parent_map <- setNames(as.integer(outbreak$parent_id), as.character(outbreak$id))
  res <- integer(0)
  cur <- id
  if (include_self) res <- c(res, cur)
  visited <- character(0)
  while (TRUE) {
    cur_chr <- as.character(cur)
    if (!cur_chr %in% names(parent_map)) break
    p <- parent_map[cur_chr]
    if (is.na(p)) break
    # cycle detection
    if (cur_chr %in% visited) {
      warning("cycle detected while walking infectors")
      break
    }
    visited <- c(visited, cur_chr)
    res <- c(res, as.integer(p))
    cur <- as.integer(p)
  }
  unique(res)
}

# Find the most recent common infector (MRCI) of two ids in the outbreak
# outbreak: data.frame with columns 'id', 'parent_id', 'infection_time', 'removal_time'
# id1, id2: integer ids
# include_self: if TRUE, include the ids themselves as the first element in their infector lists
# returns: list with elements 'mrci' (the MRCI id or NA if none) and 'distance' (the transmission distance between id1 and id2 via the M
find_mrci <- function(outbreak, id1, id2, include_self = TRUE) {
  inf1 <- find_infectors(outbreak, id1, include_self)
  inf2 <- find_infectors(outbreak, id2, include_self)
  infector <- intersect(inf1, inf2)[1]
  position1 <- which(inf1 == infector)
  position2 <- which(inf2 == infector)
  if (position1 == 1 && position2 == 1) {
    distance <- 0 # same ids
  } else if (position1 == 1) {
    infection_time <- outbreak[inf2[position2 - 1], ]$infection_time
    distance <- outbreak[id1, ]$removal_time + outbreak[id2, ]$removal_time - 2 * infection_time
  } else if (position2 == 1) {
    infection_time <- outbreak[inf1[position1 - 1], ]$infection_time
    distance <- outbreak[id1, ]$removal_time + outbreak[id2, ]$removal_time - 2 * infection_time
  } else {
    infection_time <- min(outbreak[inf1[position1 - 1], ]$infection_time, outbreak[inf2[position2 - 1], ]$infection_time)
    distance <- outbreak[id1, ]$removal_time + outbreak[id2, ]$removal_time - 2 * infection_time
  }
  list(mrci = infector, distance = ifelse(is.na(infector), NA, distance))
}

# Compute all pairwise distances using MRCA based on ancestor lists
# Returns either an n x n symmetric matrix (return_matrix=TRUE) or a long data.frame
find_all_pairwise_distances <- function(outbreak, return_matrix = TRUE) {
  if (!is.data.frame(outbreak) || !all(c("id", "parent_id", "infection_time", "removal_time") %in% names(outbreak))) {
    stop("outbreak must be a data.frame with columns 'id','parent_id','infection_time','removal_time'")
  }

  ids <- outbreak$id
  n <- length(ids)

  # maps for quick lookup
  inf_time <- setNames(as.numeric(outbreak$infection_time), as.character(outbreak$id))
  rem_time <- setNames(as.numeric(outbreak$removal_time), as.character(outbreak$id))

  # build ancestor chains for each id (include self so MRCA can be one of the ids)
  chains <- lapply(ids, function(x) find_infectors(outbreak, x, include_self = TRUE))
  names(chains) <- as.character(ids)

  # helper to compute MRCA (infector) and p_inf per pair using ancestor positions (minimize sum of positions)
  compute_mrca_and_pinf <- function(a, b) {
    ca <- chains[[as.character(a)]]
    cb <- chains[[as.character(b)]]
    if (length(ca) == 0 || length(cb) == 0) {
      return(list(mrca = NA_integer_, p_inf = NA_real_, distance = NA_real_))
    }
    common <- intersect(ca, cb)
    if (length(common) == 0) {
      return(list(mrca = NA_integer_, p_inf = NA_real_, distance = NA_real_))
    }
    # compute positions
    pos_a <- match(common, ca)
    pos_b <- match(common, cb)
    sums <- pos_a + pos_b
    best <- which.min(sums)[1]
    infector <- as.integer(common[best])
    position1 <- pos_a[best]
    position2 <- pos_b[best]

    # determine p_inf like find_mrci
    if (position1 == 1 && position2 == 1) {
      # same id
      p_inf <- inf_time[as.character(infector)]
      distance <- 0.0
    } else if (position1 == 1) {
      # infector is id a; use previous ancestor in b
      prev_b <- cb[position2 - 1]
      p_inf <- inf_time[as.character(prev_b)]
      distance <- if (is.na(p_inf)) NA_real_ else rem_time[as.character(a)] + rem_time[as.character(b)] - 2 * p_inf
    } else if (position2 == 1) {
      prev_a <- ca[position1 - 1]
      p_inf <- inf_time[as.character(prev_a)]
      distance <- if (is.na(p_inf)) NA_real_ else rem_time[as.character(a)] + rem_time[as.character(b)] - 2 * p_inf
    } else {
      prev_a <- ca[position1 - 1]
      prev_b <- cb[position2 - 1]
      p_inf <- min(inf_time[as.character(prev_a)], inf_time[as.character(prev_b)], na.rm = TRUE)
      if (is.infinite(p_inf)) p_inf <- NA_real_
      distance <- if (is.na(p_inf)) NA_real_ else rem_time[as.character(a)] + rem_time[as.character(b)] - 2 * p_inf
    }
    list(mrca = infector, p_inf = p_inf, distance = distance)
  }

  if (return_matrix) {
    mat <- matrix(NA_real_, nrow = n, ncol = n)
    rownames(mat) <- colnames(mat) <- as.character(ids)
    for (i in seq_len(n)) {
      for (j in seq.int(i, n)) {
        id1 <- ids[i]
        id2 <- ids[j]
        tmp <- compute_mrca_and_pinf(id1, id2)
        d <- tmp$distance
        mat[i, j] <- mat[j, i] <- d
      }
    }
    return(mat)
  } else {
    rows <- vector("list", n * (n - 1) / 2)
    k <- 1L
    for (i in seq_len(n)) {
      for (j in seq.int(i + 1L, n)) {
        id1 <- ids[i]
        id2 <- ids[j]
        tmp <- compute_mrca_and_pinf(id1, id2)
        rows[[k]] <- data.frame(id1 = id1, id2 = id2, common_parent = tmp$mrca, parent_infection_time = tmp$p_inf, distance = tmp$distance)
        k <- k + 1L
      }
    }
    df <- do.call(rbind, rows)
    return(df)
  }
}

# add coalescent distances to the pairwise distance matrix
add_coalescent_distance <- function(pairwise_dist, theta) {
  ncase <- dim(pairwise_dist)[1]
  coal_distance <- matrix(rexp(ncase * ncase, 1 / theta), nrow = ncase, ncol = ncase)
  coal_distance <- coal_distance + t(coal_distance)
  distance <- pairwise_dist + coal_distance
  diag(distance) <- 0.0
  distance
}

# Build a neighbor-joining tree from the outbreak distance matrix
# Requires package 'ape'
# Returns a 'phylo' object (invisible) and optionally plots the tree
build_nj_tree <- function(outbreak, root = 1, plot = FALSE, ...) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required for build_nj_tree. Please install it with install.packages('ape').")
  }
  pairwise_dist <- find_all_pairwise_distances(outbreak, return_matrix = TRUE)

  # convert to 'dist' object (rownames are ids)
  d <- stats::as.dist(pairwise_dist)
  tree <- ape::nj(d)
  tree <- ape::root(tree, root, resolve.root = TRUE)
  tree$edge.length[1] <- tree$edge.length[length(tree$edge.length)] * 0.2
  tree$edge.length[length(tree$edge.length)] <- tree$edge.length[length(tree$edge.length)] * 0.8
  if (plot) {
    ape::plot.phylo(tree, ...)
  }

  # check if the nj tree is correct
  node_dist <- cophenetic.phylo(tree)
  if (max(abs(node_dist - pairwise_dist)) > 10^-10) {
    stop("pairwise distances are not consistent with the tip distances in the nj tree")
  }

  invisible(tree)
}

# ----- Convert computeMatWIW output into a transmission edge list -----
# The object `result` returned by computeMatWIW is expected to be a matrix
# where rows are infectors (including a 0 for external/source) and columns are
# infectees (cases). Each cell holds the posterior probability that the row
# infector infected the column case. We'll take, for each column, the infector
# with the highest posterior probability (tie broken by smallest id), and
# produce an edge list (infector -> infectee).

transmission_edges_matrixwiw <- function(matrix_wiw) {
  mat <- matrix_wiw
  # rows correspond to possible infectors; ensure rownames/colnames exist
  infector_ids <- rownames(mat)
  infectee_ids <- colnames(mat)

  # If rownames/colnames are missing, create numeric ids with 0 allowed
  if (is.null(infector_ids)) infector_ids <- as.character(seq_len(nrow(mat)))
  if (is.null(infectee_ids)) infectee_ids <- as.character(seq_len(ncol(mat)))

  edges <- data.frame(infector = character(0), infectee = character(0), prob = numeric(0), stringsAsFactors = FALSE)

  for (j in seq_len(ncol(mat))) {
    col_probs <- mat[, j]
    # pick the index of the maximum posterior probability; if all NA, set NA
    if (all(is.na(col_probs))) {
      best_idx <- NA
      best_prob <- NA
    } else {
      best_idx <- which(col_probs == max(col_probs, na.rm = TRUE))
      # tie-break by smallest index
      best_idx <- best_idx[1]
      best_prob <- col_probs[best_idx]
    }

    infectee <- infectee_ids[j]
    if (is.na(best_idx) || best_prob == 0.0) {
      infector <- NA
    } else {
      infector <- infector_ids[best_idx]
    }

    edges <- rbind(edges, data.frame(infectee = infectee, infector = infector, prob = best_prob, stringsAsFactors = FALSE))
  }

  # Normalize types
  edges$infector <- as.character(edges$infector)
  edges$infectee <- as.character(edges$infectee)
  rownames(edges) <- NULL
  edges
}

sim_coaltree <- function(timetree, murate = 10^-4, theta = 10^-4) {
    # add mutations
    phytree <- timetree
    phytree$edge.length <- phytree$edge.length * murate
    phytreestring <- write.tree(phytree)

    # calculate tip-root distances
    spname <- phytree$tip.label
    nspecies <- length(spname)
    root <- nspecies + 1
    root_to_tip <- dist.nodes(phytree)[1:nspecies, root]
    diff_tip_length <- max(root_to_tip) - root_to_tip

    # convert to a clock tree by adding to tips
    clocktree <- phytree
    index <- match(spname, clocktree$tip.label)
    tip_index <- order(clocktree$edge[, 2])[index]
    clocktree$edge.length[tip_index] <- clocktree$edge.length[tip_index] + diff_tip_length
    
    # use clocktree to simulate coalescent tree
    sptree <- read.tree.nodes(write.tree(clocktree), name = spname)
    nodematrix <- sptree$nodes
    nodematrix[, 5] <- theta
    genetree <- sim.coaltree.sp(dim(nodematrix)[1], nodematrix, numcase, seq = rep(1, numcase), name = spname)$gt
    coaltree <- read.tree(text = genetree)

    # convert back to original branch lengths by subtracting from tips
    index <- match(spname, coaltree$tip.label)
    tip_index <- order(coaltree$edge[, 2])[index]    
    coaltree$edge.length[tip_index] <- coaltree$edge.length[tip_index] - diff_tip_length
    return(coaltree)
}

summary_bestnet <- function(bestnet_output_file, burnin = 0.1, onsite_time, removal_time, plot = TRUE) {
  output <- read.csv(bestnet_output_file)
  numcase <- (dim(output)[2] - 6) / 2
  transmission <- data.frame(
    id = integer(numcase),
    parent_id = integer(numcase),
    infection_time = numeric(numcase),
    latent_period = numeric(numcase),
    onsite_time = numeric(numcase),
    infectious_period = numeric(numcase),
    removal_time = numeric(numcase),
    infector_post_probability = numeric(numcase),
    stringsAsFactors = FALSE
  )
  transmission$onsite_time <- onsite_time
  transmission$removal_time <- removal_time

  # convergence plots
  if (plot) {
    par(mfrow = c(2, 2))
    plot(output$logLikelihood, type = "l", xlab = "iteration", ylab = "loglikelihood")
    plot(output$theta, type = "l", xlab = "iteration", ylab = "theta")
    plot(output$mu, type = "l", xlab = "iteration", ylab = "mutation rate")
    plot(output$infection_rate, type = "l", xlab = "iteration", ylab = "infection rate")
    par(mfrow = c(1, 1))
  }

  # parameter estimates
  post_sample <- output[(floor(dim(output)[1] * burnin) + 1):dim(output)[1], ]
  theta_est <- mean(post_sample$theta)
  mu_est <- mean(post_sample$mu)
  infrate_est <- mean(post_sample$infection_rate)

  # infector posterior probabilities
  mat <- matrix(0, numcase, numcase)
  for (i in 2:numcase) {
    x <- table(post_sample[, 6 + i])
    mat[as.numeric(names(x)), i] <- x / sum(x)
  }
  x <- transmission_edges_matrixwiw(mat)
  transmission$id <- as.numeric(x[, 1])
  transmission$parent_id <- as.numeric(x[, 2])
  transmission$infector_post_probability <- x[, 3]

  inftime <- rep(0, numcase)
  for (j in 2:numcase) {
    infection_time <- post_sample[, 6 + numcase + j]
    infector <- transmission$parent_id[j]
    index <- which(post_sample[, 6 + j] == infector)
    inftime[j] <- mean(infection_time[index])
  }
  transmission$infection_time <- inftime
  transmission$latent_period <- transmission$onsite_time - transmission$infection_time
  transmission$infectious_period <- transmission$removal_time - transmission$onsite_time
  list(theta = theta_est, mu = mu_est, infrate = infrate_est, transmission = transmission)
}

ttree_from_transmission <- function(outbreak, dateLastSample) {
  outbreak$removal_time <- outbreak$removal_time - max(outbreak$removal_time) + dateLastSample
  outbreak$onsite_time <- outbreak$removal_time - outbreak$infectious_period
  outbreak$infection_time <- outbreak$onsite_time - outbreak$latent_period

  ttree <- cbind(outbreak$infection_time, outbreak$removal_time, outbreak$parent_id)
  ttree[1, 3] <- 0
  name <- as.character(outbreak$id)
  tree <- list(ttree = ttree, nam = name)
  class(tree) <- "ttree"
  tree
}

transplot <- function(transmission, style = 1, vertex_colors = rep("lightblue", length(transmission[, 1])), vertex_sizes = rep(12, length(transmission[, 1])), vertex_label_cex = rep(1.5, length(transmission[, 1])), showlabel = TRUE, dateLastSample = 2008) {
  # Filter out zero-probability (and missing) edges for plotting only
  edges_plot <- transmission
  # treat NA as zero for plotting purposes
  edges_plot$infector_post_probability[is.na(edges_plot$infector_post_probability)] <- 0
  edges_plot <- edges_plot[edges_plot$infector_post_probability > 0, , drop = FALSE]
  if (nrow(edges_plot) == 0) {
    stop("No edges with prob > 0 to plot; skipping plotting.")
  }
  # convert NA infector (external) to "source"
  edges_plot$parent_id[is.na(edges_plot$parent_id)] <- "source"
  # create a graph
  g <- graph_from_data_frame(edges_plot[, c("parent_id", "id")], directed = TRUE)
  V(g)$name <- as.character(V(g)$name)
  target_idx <- which(V(g)$name == "1")
  if (length(target_idx) == 1) {
    vertex_colors[target_idx] <- "grey"
    vertex_sizes[target_idx] <- vertex_sizes[target_idx] + 2
    vertex_label_cex[target_idx] <- vertex_label_cex[target_idx] + 0.2
  }
  E(g)$weight <- round(edges_plot$infector_post_probability, digits = 2)
  E(g)$color <- ifelse(edges_plot$infector_post_probability < 0.5, "pink",
    ifelse(edges_plot$infector_post_probability < 0.75, "#a1a112", "#6b1a1a")
  )
  if (showlabel) {
    E(g)$label <- as.character(E(g)$weight)
  } else {
    E(g)$label <- ""
  }

  if (style == 1) { # network styple plot
    plot(
      g,
      vertex.size = vertex_sizes,
      vertex.color = vertex_colors,
      vertex.label.cex = vertex_label_cex,
      edge.label = E(g)$label, # show edge attribute
      # edge.width = E(g)$weight, # edge width proportional to weight
      # edge.label.cex = 1.5,
      edge.color = E(g)$color
    )
    legend("topright", legend = c("< 0.5", "0.5 - 0.75", "> 0.75"), col = c("pink", "#a1a112", "#6b1a1a"), pch = 19, bty = "n", cex = 1.5)
  } else if (style == 2) { # rooted-tree-style plot
    root_name <- if ("source" %in% V(g)$name) "source" else V(g)$name[1]
    layout <- layout_as_tree(g, root = which(V(g)$name == root_name))
    plot(
      g,
      layout = layout,
      vertex.size = vertex_sizes,
      vertex.color = vertex_colors,
      vertex.label.cex = vertex_label_cex,
      edge.color = E(g)$color
    )
    legend("topright", legend = c("< 0.5", "0.5 - 0.75", "> 0.75"), col = c("pink", "#a1a112", "#6b1a1a"), pch = 19, bty = "n", cex = 1.5)
  } else if (style == 3) { # timeline style plot
    ttree <- ttree_from_transmission(transmission, dateLastSample = dateLastSample)
    plot(ttree)
  } else if (style == 4) { # transmission style plot
    ttree <- ttree_from_transmission(transmission, dateLastSample = dateLastSample)
    plot(ttree, type = "detailed", w.shape = 10, w.scale = 0.1)
  }
}
