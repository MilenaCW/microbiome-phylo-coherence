consentrait_signed <- function(
  tree,
  trait_values,
  frac_consensus       = 0.9,
  n_shuffles           = 100,
  singleton_depth_frac = 0,
  weight_clades        = FALSE,
  seed                 = 1
) {
  
  if (!is.null(seed)) set.seed(seed)

  # ---- input validation ----
  if (!ape::is.rooted(tree))
    warning("tree is not rooted; results may be unreliable")
  tv_unique <- unique(trait_values[!is.na(trait_values)])
  if (!all(tv_unique %in% c(-1, 1)))
    stop("trait_values must contain only +1 or -1")
  if (is.null(tree$edge.length))
    tree$edge.length <- rep(1, nrow(tree$edge))

  # ---- tree preprocessing ----
  tree <- ape::multi2di(tree)
  tree$edge.length[is.na(tree$edge.length) | tree$edge.length <= 0] <- 1e-6

  # ---- align tip labels ----
  common       <- intersect(tree$tip.label, names(trait_values))
  if (length(common) == 0) stop("No tips in common between tree and trait_values")
  n_drop_tree  <- length(tree$tip.label) - length(common)
  n_drop_trait <- length(trait_values)   - length(common)
  if (n_drop_tree  > 0) {
    warning(n_drop_tree, " tips in tree not in trait_values; dropped from tree")
    tree <- ape::keep.tip(tree, common)
  }
  if (n_drop_trait > 0)
    warning(n_drop_trait, " names in trait_values not in tree; ignored")
  trait_values <- trait_values[tree$tip.label]  # aligned to tree tip order

  n_tips  <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  root_nd <- n_tips + 1L

  # ---- pre-order pass: distance of each node from root ----
  tree_cl  <- ape::reorder.phylo(tree, "cladewise")
  cl_p     <- tree_cl$edge[, 1L]
  cl_ch    <- tree_cl$edge[, 2L]
  cl_len   <- tree_cl$edge.length
  nd_depth <- numeric(n_nodes)  # root depth = 0

  for (i in seq_along(cl_p))
    nd_depth[cl_ch[i]] <- nd_depth[cl_p[i]] + cl_len[i]

  # ---- post-order pass: subtree tip counts and sum of tip depths ----
  # sub_tds[node] = sum of nd_depth[tip] for all tips in node's subtree
  # sub_n[node]   = number of tips in node's subtree
  tree_po <- ape::reorder.phylo(tree, "postorder")
  po_p  <- tree_po$edge[, 1L]
  po_ch <- tree_po$edge[, 2L]

  sub_n   <- integer(n_nodes)
  sub_tds <- numeric(n_nodes)
  sub_n  [1:n_tips] <- 1L
  sub_tds[1:n_tips] <- nd_depth[1:n_tips]

  for (i in seq_along(po_p)) {
    sub_n  [po_p[i]] <- sub_n  [po_p[i]] + sub_n  [po_ch[i]]
    sub_tds[po_p[i]] <- sub_tds[po_p[i]] + sub_tds[po_ch[i]]
  }

  # terminal branch lengths, indexed by tip number (for singletons)
  tip_edge <- match(seq_len(n_tips), tree$edge[, 2L])
  term_bl  <- tree$edge.length[tip_edge]

  # ---- helper: post-order sub_sum from a (possibly shuffled) trait vector ----
  # sub_sum[node] = sum of trait values for all tips in node's subtree
  make_sub_sum <- function(vals) {
    ss           <- numeric(n_nodes)
    ss[1:n_tips] <- vals
    for (i in seq_along(po_p))
      ss[po_p[i]] <- ss[po_p[i]] + ss[po_ch[i]]
    ss
  }

  # ---- helper: pre-order coherence pass ----
  # Coherence criterion: abs(sub_sum[node]) / sub_n[node] > 2*frac_consensus - 1
  # Equivalent to: max(n_plus, n_minus) / n_tips_in_clade > frac_consensus
  thresh <- 2 * frac_consensus - 1

  run_pass <- function(sub_sum_vec) {
    skip <- logical(n_nodes)  # TRUE = inside an already-claimed coherent clade

    # Pre-allocate output vectors: worst case is n_tips records (all singletons)
    rec_dir  <- numeric(n_tips)
    rec_dep  <- numeric(n_tips)
    rec_sz   <- integer(n_tips)
    rec_type <- character(n_tips)
    rec_node <- integer(n_tips)
    n_rec    <- 0L

    add_rec <- function(dir, dep, sz, type, nd) {
      n_rec <<- n_rec + 1L
      rec_dir [n_rec] <<- dir   # direction: +1 or -1 for clades; actual trait value for singletons
      rec_dep [n_rec] <<- dep   # depth: mean tip depth for clades (depth of tips to root - depth of node to root); terminal branch length for singletons
      rec_sz  [n_rec] <<- sz    # size: number of tips in clade; always 1 for singletons
      rec_type[n_rec] <<- type  # type: "clade" or "singleton"
      rec_node[n_rec] <<- nd    # node number in tree: internal node for clades; tip number for singletons
    }

    # Check root separately (root never appears as a child in tree$edge)
    if (abs(sub_sum_vec[root_nd]) / sub_n[root_nd] > thresh) {
      add_rec(sign(sub_sum_vec[root_nd]),
              sub_tds[root_nd] / sub_n[root_nd],  # nd_depth[root] = 0
              sub_n[root_nd], "clade", root_nd)
      skip[root_nd] <- TRUE
    }

    # Pre-order edge traversal: parents evaluated before children
    for (i in seq_along(cl_p)) {
      p  <- cl_p[i]; ch <- cl_ch[i]
      # Propagate skip from parent to child
      if (skip[p]) { skip[ch] <- TRUE; next }
      if (ch <= n_tips) next  # tip nodes: handled in singleton pass below

      if (abs(sub_sum_vec[ch]) / sub_n[ch] > thresh) {
        add_rec(sign(sub_sum_vec[ch]),
                sub_tds[ch] / sub_n[ch] - nd_depth[ch],
                sub_n[ch], "clade", ch)
        skip[ch] <- TRUE
      }
    }

    # Singletons: tips not inside any coherent clade
    # depth = terminal branch length; weight_singletons controls contribution to tau_D
    for (tip in seq_len(n_tips)) {
      if (!skip[tip])
        add_rec(sub_sum_vec[tip], term_bl[tip], 1L, "singleton", tip)
    }

    # Trim pre-allocated vectors to actual length and build data frame
    if (n_rec == 0L) {
      df <- data.frame(
        direction = numeric(0), depth = numeric(0), size = integer(0),
        type = character(0), node = integer(0), stringsAsFactors = FALSE
      )
    } else {
      df <- data.frame(
        direction = rec_dir [1:n_rec],
        depth     = rec_dep [1:n_rec],
        size      = rec_sz  [1:n_rec],
        type      = rec_type[1:n_rec],
        node      = rec_node[1:n_rec],
        stringsAsFactors = FALSE
      )
    }

    # Compute tau_D: mean depth across clades (and optionally singletons).
    # singleton_depth_frac scales the singleton terminal branch before averaging;
    # each singleton still counts as one entry (does not change the denominator).
    # weight_clades makes clade size the weight in a true weighted mean.
    is_cl <- df$type == "clade"
    is_sg <- df$type == "singleton"
    cw    <- if (weight_clades) df$size[is_cl] else rep(1L, sum(is_cl))
    if (singleton_depth_frac > 0 && any(is_sg)) {
      all_d <- c(df$depth[is_cl], singleton_depth_frac * df$depth[is_sg])
      all_w <- c(cw, rep(1L, sum(is_sg)))
    } else {
      all_d <- df$depth[is_cl]
      all_w <- cw
    }
    tau_D <- if (sum(all_w) == 0) NA_real_ else stats::weighted.mean(all_d, all_w)

    list(clades = df, tau_D = tau_D)
  }

  # ---- true run ----
  true_res <- run_pass(make_sub_sum(trait_values))

  # ---- null runs: shuffle tip labels, recompute sub_sum only ----
  null_tau_D <- numeric(n_shuffles)
  null_recs  <- vector("list", n_shuffles)

  for (k in seq_len(n_shuffles)) {
    nr            <- run_pass(make_sub_sum(sample(trait_values)))
    null_tau_D[k] <- nr$tau_D
    if (nrow(nr$clades) > 0) {
      nr$clades$shuffle_id <- k
      null_recs[[k]] <- nr$clades
    }
  }

  null_clades <- do.call(rbind, Filter(Negate(is.null), null_recs))
  if (is.null(null_clades))
    null_clades <- data.frame(
      direction = numeric(0), depth = numeric(0), size = integer(0),
      type = character(0), node = integer(0), shuffle_id = integer(0),
      stringsAsFactors = FALSE
    )

  list(
    tau_D       = true_res$tau_D,
    clades      = true_res$clades,
    null_tau_D  = null_tau_D,
    null_clades = null_clades
  )
}
