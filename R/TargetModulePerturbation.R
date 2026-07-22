# Targeted extension of ModulePerturbation().

.mpt_log_normalize <- function(X, col_sums, scale.factor = 1e4) {
  if (any(!is.finite(col_sums)) || any(col_sums <= 0)) {
    stop("Every cell must have a positive, finite library size.")
  }

  if (inherits(X, "sparseMatrix")) {
    X <- methods::as(X, "CsparseMatrix")
    inv_col_sums <- scale.factor / col_sums
    X@x <- X@x * rep(inv_col_sums, diff(X@p))
    X@x <- log1p(X@x)
    return(X)
  }

  log1p(sweep(X, 2, col_sums, "/") * scale.factor)
}

.mpt_get_assay_data <- function(seurat_obj, assay, layer, slot) {
  if (hdWGCNA::CheckSeurat5()) {
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  } else {
    Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  }
}

.mpt_set_data_layer <- function(seurat_obj, assay, value) {
  if (hdWGCNA::CheckSeurat5()) {
    Seurat::SetAssayData(
      seurat_obj,
      assay = assay,
      layer = "data",
      new.data = value
    )
  } else {
    Seurat::SetAssayData(
      seurat_obj,
      assay = assay,
      slot = "data",
      new.data = value
    )
  }
}

.mpt_targeted_transitions <- function(
    seurat_obj,
    perturbation_name,
    features,
    graph,
    source_cells,
    corr_sigma = 0.05,
    assay = "RNA"
) {
  if (!is.numeric(corr_sigma) || length(corr_sigma) != 1L ||
      !is.finite(corr_sigma) || corr_sigma <= 0) {
    stop("corr_sigma must be one positive, finite number.")
  }

  all_cells <- colnames(seurat_obj)
  source_cells <- intersect(all_cells, source_cells)
  non_source_cells <- setdiff(all_cells, source_cells)

  # Rows are possible destinations and columns are transition sources.
  # Keep every row, but remove all outgoing edges from non-targeted cells.
  cell_graph <- Matrix::Matrix(seurat_obj@graphs[[graph]], sparse = TRUE)
  cell_graph <- cell_graph[all_cells, all_cells]
  diag(cell_graph) <- 1
  if (length(non_source_cells) > 0L) {
    cell_graph[, non_source_cells] <- 0
  }
  cell_graph <- methods::as(cell_graph, "CsparseMatrix")

  exp_obs <- .mpt_get_assay_data(
    seurat_obj, assay = assay, layer = "data", slot = "data"
  )
  exp_per <- .mpt_get_assay_data(
    seurat_obj, assay = perturbation_name, layer = "data", slot = "data"
  )

  features <- intersect(features, intersect(rownames(exp_obs), rownames(exp_per)))
  if (length(features) < 2L) {
    stop("At least two shared transition features are required.")
  }

  exp_obs <- methods::as(exp_obs[features, all_cells, drop = FALSE], "CsparseMatrix")
  delta <- methods::as(
    exp_per[features, all_cells, drop = FALSE] - exp_obs,
    "CsparseMatrix"
  )

  cc <- SparseColDeltaCor(exp_obs, delta, cell_graph)
  diag(cc) <- 0

  tp <- exp(cc / corr_sigma) * cell_graph
  col_sums <- Matrix::colSums(tp)

  # A valid targeted source has a self-edge and should therefore have a
  # positive column sum. Guard against numerical pathologies nevertheless.
  bad_source <- source_cells[
    !is.finite(col_sums[source_cells]) | col_sums[source_cells] <= 0
  ]
  if (length(bad_source) > 0L) {
    warning(
      length(bad_source),
      " targeted source cell(s) had an invalid transition column; ",
      "their transition probabilities were set to zero.",
      call. = FALSE
    )
    tp[, bad_source] <- 0
    col_sums[bad_source] <- 1
  }

  if (length(source_cells) > 0L) {
    tp[, source_cells] <- t(
      t(tp[, source_cells, drop = FALSE]) / col_sums[source_cells]
    )
  }
  if (length(non_source_cells) > 0L) {
    tp[, non_source_cells] <- 0
  }

  tp <- methods::as(tp, "dgCMatrix")
  rownames(tp) <- all_cells
  colnames(tp) <- all_cells
  seurat_obj@graphs[[paste0(perturbation_name, "_tp")]] <- Seurat::as.Graph(tp)

  seurat_obj
}

#' Perturb one co-expression module in a selected population only
#'
#' A targeted counterpart to \code{ModulePerturbation}. Direct
#' perturbation and gene-network propagation are restricted to cells selected
#' by \code{target.by} and \code{target_name}. All other cells remain expression-identical to
#' baseline. Cell-transition probabilities use the original shared graph, with
#' targeted cells as sources and all cells retained as possible destinations.
#'
#' @param seurat_obj Seurat object containing an hdWGCNA experiment.
#' @param mod Co-expression module to perturb.
#' @param perturb_dir Perturbation direction/magnitude passed to
#'   \code{ApplyPerturbation}.
#' @param perturbation_name Name of the new perturbation assay.
#' @param target.by Metadata column defining the cells eligible for perturbation.
#' @param target_name One or more values in target.by to perturb.
#' @param graph Existing shared cell-cell graph, for example "RNA_snn".
#' @param group.by Optional metadata column used only for ZINB model fitting
#'   within the targeted cells. It does not select cells.
#' @param n_hubs Number of module hub genes receiving the primary perturbation.
#' @param perturb_mode "zinb" or "multiplicative".
#' @param n_iters Number of gene-network propagation iterations.
#' @param expand_module Number of high-kME grey genes added to a small module.
#' @param delta_scale Per-iteration propagation scale.
#' @param row_normalize Whether to row-normalize the gene network.
#' @param prune_network Whether to prune weak gene-network edges.
#' @param prune_percentile Edge-weight percentile used when pruning.
#' @param corr_sigma Transition correlation softmax scale.
#' @param use_velocyto Retained for call compatibility. Targeted transitions
#'   currently require FALSE so the shared graph can be source-masked before
#'   the sparse correlation calculation.
#' @param layer Baseline expression layer used for primary perturbation.
#' @param slot Seurat-v4 equivalent of layer.
#' @param assay Baseline assay.
#' @param n_workers Workers used by the ZINB primary perturbation.
#' @param custom_network Optional gene-by-gene network.
#' @param custom_modules Optional data.frame with gene_name and module columns.
#' @param custom_weights Optional module-table column used to rank hub genes.
#' @param wgcna_name hdWGCNA experiment name; defaults to active_wgcna.
#'
#' @return The input Seurat object with a targeted perturbation assay, a
#'   \code{<perturbation_name>_tp} transition graph, a logical source-cell metadata
#'   column, and targeting provenance in \code{@misc$targeted_perturbations}.
#'
#' @seealso \code{\link{ModulePerturbation}}, \code{\link{ApplyPerturbation}},
#'   \code{\link{ApplyPropagation}}
#' @import Seurat
#' @export
TargetModulePerturbation <- function(
    seurat_obj,
    mod,
    perturb_dir,
    perturbation_name,
    target.by,
    target_name,
    graph = "RNA_snn",
    group.by = NULL,
    n_hubs = 5,
    perturb_mode = "zinb",
    n_iters = 3,
    expand_module = 0,
    delta_scale = 0.2,
    row_normalize = FALSE,
    prune_network = FALSE,
    prune_percentile = 0.95,
    corr_sigma = 0.05,
    use_velocyto = FALSE,
    layer = "counts",
    slot = "counts",
    assay = "RNA",
    n_workers = 1,
    custom_network = NULL,
    custom_modules = NULL,
    custom_weights = NULL,
    wgcna_name = NULL
) {
  if (!requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("hdWGCNA", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE)) {
    stop("Seurat, hdWGCNA, and Matrix must be installed and loadable.")
  }

  perturbation_name <- .sanitize_perturbation_name(perturbation_name)
  perturb_mode <- match.arg(perturb_mode, c("zinb", "multiplicative"))

  if (!identical(use_velocyto, FALSE)) {
    stop(
      "TargetModulePerturbation currently requires use_velocyto = FALSE ",
      "so transition-source cells can be masked safely."
    )
  }

  if (!(assay %in% names(seurat_obj@assays))) {
    stop("Assay '", assay, "' was not found.")
  }
  if (perturbation_name %in% names(seurat_obj@assays)) {
    stop("Assay '", perturbation_name, "' already exists; choose a new name.")
  }
  if (!(graph %in% names(seurat_obj@graphs))) {
    stop("Graph '", graph, "' was not found.")
  }
  if (!(target.by %in% colnames(seurat_obj@meta.data))) {
    stop("target.by column '", target.by, "' was not found in metadata.")
  }
  if (!is.null(group.by) && !(group.by %in% colnames(seurat_obj@meta.data))) {
    stop("group.by column '", group.by, "' was not found in metadata.")
  }
  if (!is.numeric(perturb_dir) || length(perturb_dir) != 1L ||
      !is.finite(perturb_dir)) {
    stop("perturb_dir must be one finite numeric value.")
  }
  if (!is.numeric(n_hubs) || length(n_hubs) != 1L || n_hubs < 1) {
    stop("n_hubs must be a positive integer.")
  }

  target_mask <- !is.na(seurat_obj@meta.data[[target.by]]) &
    seurat_obj@meta.data[[target.by]] %in% target_name
  target_cells <- rownames(seurat_obj@meta.data)[target_mask]
  non_target_cells <- setdiff(colnames(seurat_obj), target_cells)
  if (length(target_cells) == 0L) {
    stop("No cells matched target.by = '", target.by, "' and target_name.")
  }
  if (length(non_target_cells) == 0L) {
    warning("Every cell is targeted; no unchanged reference population remains.")
  }

  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  if (!is.null(custom_network) || !is.null(custom_modules)) {
    if (is.null(custom_network) || is.null(custom_modules)) {
      stop("custom_network and custom_modules must be supplied together.")
    }
    if (!all(c("gene_name", "module") %in% colnames(custom_modules))) {
      stop("custom_modules must contain gene_name and module columns.")
    }
    message("Using provided custom network and module table...")
    valid_genes <- intersect(rownames(custom_network), custom_modules$gene_name)
    net <- custom_network[valid_genes, valid_genes, drop = FALSE]
    modules <- custom_modules[custom_modules$gene_name %in% valid_genes, , drop = FALSE]
  } else {
    net <- hdWGCNA::GetTOM(seurat_obj, wgcna_name)
    modules <- hdWGCNA::GetModules(seurat_obj, wgcna_name)
  }

  cur_mod <- modules[modules$module == mod, , drop = FALSE]
  if (nrow(cur_mod) == 0L) {
    stop("Module '", mod, "' was not found.")
  }
  module_genes <- unique(as.character(cur_mod$gene_name))

  if (!is.null(custom_weights)) {
    if (!(custom_weights %in% colnames(cur_mod))) {
      stop("custom_weights column '", custom_weights, "' was not found.")
    }
    message("Ranking hubs using custom weight column: ", custom_weights)
    cur_mod <- cur_mod[order(cur_mod[[custom_weights]], decreasing = TRUE), , drop = FALSE]
    hub_genes <- head(as.character(cur_mod$gene_name), n_hubs)
  } else if (any(grepl("kME", colnames(cur_mod)))) {
    message("Ranking hubs using kME...")
    hubs <- hdWGCNA::GetHubGenes(
      seurat_obj, n = n_hubs, wgcna_name = wgcna_name
    )
    hub_genes <- as.character(hubs$gene_name[hubs$module == mod])

    if (expand_module > 0L) {
      kme_col <- paste0("kME_", mod)
      if (kme_col %in% colnames(modules)) {
        grey <- modules[modules$module == "grey", , drop = FALSE]
        grey <- grey[order(grey[[kme_col]], decreasing = TRUE), , drop = FALSE]
        module_genes <- unique(c(
          module_genes,
          head(as.character(grey$gene_name), expand_module)
        ))
      }
    }
  } else {
    message("No weight column found. Calculating hub genes based on network degree (column sums)...")
    mod_net <- net[module_genes, module_genes, drop = FALSE]
    degree <- Matrix::colSums(mod_net)
    hub_genes <- names(sort(degree, decreasing = TRUE))[seq_len(
      min(n_hubs, length(degree))
    )]
  }

  exp_obs <- .mpt_get_assay_data(
    seurat_obj, assay = assay, layer = layer, slot = slot
  )
  module_genes <- intersect(module_genes, intersect(rownames(net), rownames(exp_obs)))
  hub_genes <- intersect(hub_genes, module_genes)
  if (length(hub_genes) == 0L) {
    stop("No selected hub genes were found in both the network and expression matrix.")
  }
  if (length(module_genes) < 2L) {
    stop("At least two module genes are required for propagation/transitions.")
  }

  message(
    "Targeting ", length(target_cells), " cell(s); preserving ",
    length(non_target_cells), " cell(s) unchanged."
  )
  message("Selected Hubs: ", paste(hub_genes, collapse = ", "))

  message("Applying primary in-silico perturbation to hub genes...")
  # ApplyPerturbation's ZINB grouping currently assumes that all metadata cells
  # are present in its expression matrix. Supplying a target-only Seurat object
  # makes that assumption true without changing ApplyPerturbation().
  target_obj <- subset(seurat_obj, cells = target_cells)
  Seurat::DefaultAssay(target_obj) <- assay
  target_exp <- exp_obs[, target_cells, drop = FALSE]
  target_exp_per <- ApplyPerturbation(
    target_obj,
    target_exp,
    features = hub_genes,
    perturb_dir = perturb_dir,
    perturb_mode = perturb_mode,
    cells_use = target_cells,
    group.by = group.by,
    layer = layer,
    slot = slot,
    assay = assay,
    n_workers = n_workers
  )

  exp_per <- exp_obs
  exp_per[, target_cells] <- target_exp_per[, target_cells, drop = FALSE]
  exp_per <- methods::as(exp_per, "CsparseMatrix")

  library_sizes <- Matrix::colSums(exp_obs)
  log_obs <- .mpt_log_normalize(exp_obs, library_sizes)
  log_primary <- .mpt_log_normalize(exp_per, library_sizes)
  delta_log <- methods::as(
    log_primary[module_genes, , drop = FALSE] -
      log_obs[module_genes, , drop = FALSE],
    "CsparseMatrix"
  )

  cur_net <- net[module_genes, module_genes, drop = FALSE]
  message("Applying log-space signal propagation throughout co-expression network...")
  log_simulated_mod <- ApplyPropagation(
    log_obs_mod = log_obs[module_genes, , drop = FALSE],
    delta_log = delta_log,
    network = cur_net,
    n_iters = n_iters,
    delta_scale = delta_scale,
    row_normalize = row_normalize,
    prune_network = prune_network,
    prune_percentile = prune_percentile
  )

  non_module_genes <- setdiff(rownames(exp_obs), module_genes)
  log_simulated <- rbind(
    log_obs[non_module_genes, , drop = FALSE],
    methods::as(log_simulated_mod, "CsparseMatrix")
  )[rownames(exp_obs), colnames(exp_obs), drop = FALSE]

  # Enforce the scientific invariant explicitly, rather than relying only on
  # the fact that zero non-target deltas remain zero under gene multiplication.
  if (length(non_target_cells) > 0L) {
    exp_per[, non_target_cells] <- exp_obs[, non_target_cells, drop = FALSE]
    log_simulated[, non_target_cells] <- log_obs[, non_target_cells, drop = FALSE]
  }

  seurat_obj[[perturbation_name]] <- Seurat::CreateAssayObject(counts = exp_per)
  seurat_obj <- .mpt_set_data_layer(
    seurat_obj, assay = perturbation_name, value = log_simulated
  )

  message("Computing cell-cell transition probabilities based on the perturbation...")
  seurat_obj <- .mpt_targeted_transitions(
    seurat_obj,
    perturbation_name = perturbation_name,
    features = module_genes,
    graph = graph,
    source_cells = target_cells,
    corr_sigma = corr_sigma,
    assay = assay
  )

  source_col <- make.names(paste0(perturbation_name, "_source"))
  seurat_obj@meta.data[[source_col]] <- rownames(seurat_obj@meta.data) %in% target_cells

  if (is.null(seurat_obj@misc$targeted_perturbations)) {
    seurat_obj@misc$targeted_perturbations <- list()
  }
  seurat_obj@misc$targeted_perturbations[[perturbation_name]] <- list(
    module = mod,
    hub_genes = hub_genes,
    module_genes = module_genes,
    target.by = target.by,
    target_name = target_name,
    source_metadata_column = source_col,
    n_target_cells = length(target_cells),
    n_unchanged_cells = length(non_target_cells),
    graph = graph,
    transition_graph = paste0(perturbation_name, "_tp"),
    params = list(
      perturb_dir = perturb_dir,
      perturb_mode = perturb_mode,
      n_hubs = n_hubs,
      n_iters = n_iters,
      delta_scale = delta_scale,
      row_normalize = row_normalize,
      prune_network = prune_network,
      prune_percentile = prune_percentile,
      corr_sigma = corr_sigma,
      group.by = group.by,
      assay = assay,
      layer = layer
    )
  )

  # Final hard checks: all non-target counts and normalized values must be
  # exactly identical to baseline, and non-target transition columns must be 0.
  if (length(non_target_cells) > 0L) {
    if (Matrix::nnzero(exp_per[, non_target_cells] - exp_obs[, non_target_cells]) != 0L) {
      stop("Internal invariant failed: non-target counts changed.")
    }
    if (Matrix::nnzero(log_simulated[, non_target_cells] - log_obs[, non_target_cells]) != 0L) {
      stop("Internal invariant failed: non-target normalized expression changed.")
    }
    tp <- seurat_obj@graphs[[paste0(perturbation_name, "_tp")]]
    if (Matrix::nnzero(tp[, non_target_cells]) != 0L) {
      stop("Internal invariant failed: non-target transition-source columns are nonzero.")
    }
  }

  seurat_obj
}

# Example for the discussed DMSO/Bortezomib experiment:
#
# seurat_sub <- TargetModulePerturbation(
#   seurat_sub,
#   mod = cur_mod,
#   perturb_dir = 1,
#   perturbation_name = make.names(paste0(cur_mod, "_up_DMSO")),
#   target.by = "drug",
#   target_name = "DMSO_TF",
#   graph = "RNA_snn",
#   n_hubs = 10,
#   delta_scale = 1.0,
#   row_normalize = FALSE,
#   use_velocyto = FALSE,
#   n_workers = 1
# )
