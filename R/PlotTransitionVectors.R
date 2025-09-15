
#' Plot Transition Vectors on a Reduced Dimensional Embedding
#'
#' This function visualizes cell-state transitions in a Seurat object by plotting transition vectors on a 
#' dimensionality-reduced embedding (e.g., UMAP). It computes transition vectors based on a specified 
#' perturbation and overlays them on a scatter plot of the cells, allowing insights into cell-state shifts
#' following the perturbation.
#'
#' @param seurat_obj A Seurat object containing single-cell data, including metadata and embeddings.
#' @param perturbation_name A character string specifying the name of the perturbation. This is used 
#'   to calculate transition vectors based on the perturbed assay.
#' @param color.by A character string indicating the metadata column to color the scatter plot points by.
#' @param reduction A character string specifying the name of the dimensional reduction to use for the 
#'   embedding (default: 'umap').
#' @param n_threads An integer specifying the number of threads to use for parallel processing 
#'   (default: 4).
#' @param arrow_scale A numeric value to scale the transition vectors, adjusting their length on the plot 
#'   (default: 1).
#' @param grid_resolution An integer specifying the resolution of the grid for overlaying transition 
#'   vectors. A higher resolution yields a finer grid, which can provide more detailed vector placement 
#'   (default: 25).
#' @param max_pct A numeric value between 0 and 1 indicating the maximum percentile of vector lengths 
#'   to display. This parameter helps limit the display of outlier vector lengths, making the plot 
#'   easier to interpret (default: 0.90).
#' @param point_alpha A numeric value controlling the transparency of scatter plot points, with 1 being 
#'   fully opaque and 0 being fully transparent (default: 0.25).
#' @param point_size A numeric value defining the size of the scatter plot points (default: 1).
#' 
#' @return A ggplot object showing the scatter plot of cells colored by the specified metadata and 
#'   overlaid with transition vectors indicating cell-state shifts.
#' 
#' @details This function retrieves the embedding coordinates from the specified dimensional reduction 
#'   in the Seurat object. It then calculates the transition vectors using `PerturbationVectors`, which 
#'   returns vector coordinates (`ars`) and distances (`arsd`). A grid of arrows is created using 
#'   `grid_vectors`, allowing for visual simplification of transitions by grouping vectors into a grid 
#'   at the specified resolution. The resulting plot includes the scatter plot of cells with color 
#'   defined by `color.by` and arrows showing the direction and relative magnitude of transitions.
#'
#' @export
PlotTransitionVectors <- function(
    seurat_obj, 
    perturbation_name,
    color.by,
    reduction = 'umap',
    n_threads = 4,
    arrow_scale = 1,
    grid_resolution = 25,
    max_pct = 0.90,
    point_alpha = 0.25,
    point_size = 1,
    arrow_size = 0.25,
    raster_dpi = 300,
    arrow_alpha = TRUE
){

    # This function gives us the dataframe with vector coordinates 
    # (ars) and vector distances (arsd)
    vectors <- PerturbationVectors(
        seurat_obj,
        perturbation_name = perturbation_name,
        reduction = reduction, 
        arrow_scale = arrow_scale,
        max_pct = max_pct
    )
    ars <- vectors$ars 
    arsd <- vectors$arsd

    # Get the UMAP embedding from the seurat object 
    # and subset to only keep cells that we have arrows for.
    emb <- Reductions(seurat_obj, reduction)@cell.embeddings[,1:2]
    colnames(emb) <- paste0('emb_', 1:ncol(emb))

    # make a grid of arrows 
    grid.df <- grid_vectors(emb[rownames(ars),], arsd, resolution=grid_resolution)

    # calculate the length of the arrow
    grid.df$vector.length <- sqrt((grid.df$end.xd - grid.df$start.emb_1)^2 + (grid.df$end.yd - grid.df$start.emb_2)^2)

    # make a dataframe to plot with ggplot
    plot_df <- as.data.frame(emb)
    plot_df$value <- seurat_obj@meta.data[,color.by]

    # arrange by value:
    plot_df <- plot_df %>% arrange(value)

    # plot the scatter plot
    p <- ggplot(plot_df, aes(x=emb_1, y=emb_2, color=value)) +
        ggrastr::rasterise(
            geom_point(alpha=point_alpha, size=point_size),
            dpi = raster_dpi
        )

    if(arrow_alpha){
        p <- p +
        geom_segment(
            data = grid.df,
            inherit.aes=FALSE,
            aes(
                x=start.emb_1, 
                y=start.emb_2, 
                xend=end.xd, 
                yend=end.yd,
                alpha = vector.length 
            ), 
            arrow=grid::arrow(length=unit(0.1, "cm")), 
            size=arrow_size
        )
    } else{
        p <- p +
        geom_segment(
            data = grid.df,
            inherit.aes=FALSE,
            aes(
                x=start.emb_1, 
                y=start.emb_2, 
                xend=end.xd, 
                yend=end.yd
            ), 
            arrow=grid::arrow(length=unit(0.1, "cm")), 
            size=arrow_size
        )
    }

    # theme elements:
    p <- p + 
        hdWGCNA::umap_theme()
    
    p

}


#' Calculate Perturbation Transition Vectors for Cell Embeddings
#'
#' This function calculates transition vectors for cells in a Seurat object based on a specified 
#' perturbation, providing insight into cell-state shifts on a dimensionality-reduced embedding.
#'
#' @param seurat_obj A Seurat object containing single-cell data, including embeddings and graphs.
#' @param perturbation_name A character string specifying the name of the perturbation. This is used to 
#'   retrieve the perturbation-specific transition probability graph.
#' @param reduction A character string specifying the name of the dimensional reduction to use for the 
#'   embedding (default: 'umap').
#' @param n_threads An integer specifying the number of threads to use for parallel processing 
#'   (default: 4).
#' @param arrow_scale A numeric value to scale the transition vectors, adjusting their length on the plot 
#'   (default: 1).
#' @param max_pct A numeric value between 0 and 1 indicating the maximum percentile of vector lengths 
#'   to display. This parameter caps the extreme values, limiting outlier vector lengths for better 
#'   visualization (default: 0.90).
#' 
#' @return A list containing two data frames:
#'   \item{ars}{A data frame with coordinates for plotting transition vectors, including initial points 
#'     (`x0`, `y0`) and final points (`x1`, `y1`) for each cell.}
#'   \item{arsd}{A data frame with the transition vector distances (`xd`, `yd`) for each cell, 
#'     normalized based on `max_pct`.}
#' 
#' @details This function retrieves the 2D embedding coordinates for the specified reduction in the Seurat 
#'   object. It then extracts the perturbation-specific transition probability graph from the Seurat object 
#'   and uses it to calculate transition vectors with the `embArrows_velocyto` helper function. Vector 
#'   lengths are scaled by `arrow_scale` and capped at the specified `max_pct` to handle extreme values.
#'   Cells with `NA` values in either `ars` or `arsd` are excluded from the final output.
#'
#' @export
PerturbationVectors <- function(
    seurat_obj,
    perturbation_name,
    reduction = 'umap',
    n_threads=4,
    arrow_scale = 1,
    max_pct = 0.90
){

    # TODO: check reduction 
    graph_name <- paste0(perturbation_name, '_tp')

    # get the 2D embedding
    emb <- Reductions(seurat_obj, reduction)@cell.embeddings[,1:2]

    # get the graph from the seurat object
    tp <- Graphs(seurat_obj, graph_name)

    # run the helper function
    arsd <- data.frame(t(embArrows_velocyto(
        emb,
        tp,
        arrow_scale,
        n_threads
    )))

    # replace NA with 0 
    arsd[is.na(arsd)] <- 0

    # format the dataframes
    rownames(arsd) <- rownames(emb)
    ars <- data.frame(cbind(emb,emb+arsd))
    colnames(ars) <- c('x0','y0','x1','y1')
    colnames(arsd) <- c('xd','yd')
    rownames(ars) <- rownames(emb)
        
    # exclude NAs:
   #  ars <- na.omit(ars)
    # arsd <- na.omit(arsd)

    max_vect <- quantile(abs(as.matrix(arsd)), max_pct)
    arsd[arsd > max_vect] <- max_vect
    arsd[arsd < -max_vect] <- -max_vect
    arsd <- arsd / max_vect  

    # return 
    list(ars=ars, arsd=arsd)

}


#' @export
#' @importFrom S4Vectors selfmatch DataFrame
#' @importFrom stats median
grid_vectors <- function(x, embedded, resolution=40, scale=TRUE, group_min=5) {
    limits <- apply(x, 2, range)
    intercept <- limits[1,]
    delta <- (limits[2,] - intercept)/resolution

    categories <- t((t(x) - intercept)/delta)
    storage.mode(categories) <- "integer"
    categories <- DataFrame(categories)
    grp <- selfmatch(categories)

    tab <- as.integer(table(grp))
    pos <- rowsum(x, grp)/tab
    vec <- rowsum(embedded, grp)/tab

    if (scale) {
        d <- sqrt(rowSums(vec^2))
        target <- sqrt(sum(delta^2))
        vec <- vec * target/median(d)/2 # divide by 2 to avoid running over into another block.
    }

    out <- data.frame(
        start=pos, 
        end=pos + vec
    )
    
    # remove outliers
    rows_keep <- which(tab >= group_min)
    out <- out[rows_keep,]
    out
    
}




embArrows <- function(embedding, transition_matrix, knn_graph, scale = 1, threshold = 0.1) {
        
    # Number of cells and dimensions
    n_cells <- nrow(embedding)
    n_dims <- ncol(embedding)
  
    # Example sparse matrix for transition probabilities
    tpb <- Matrix::Matrix(data = knn_graph, sparse = TRUE)

    # Column sums for non-zero elements
    col_sums <- colSums(tpb)

    tpb <- tpb %*% Diagonal(x = 1 / col_sums)

    # Initialize delta_embedding matrix
    delta_embedding <- matrix(0, nrow = n_cells, ncol = n_dims)
    
    # Transpose embedding for easier manipulation
    transposed_embedding <- t(embedding)
    
    # Zero vector for normalization
    zero_vec <- rep(0, n_dims)
    
    # Loop over each cell
    for (i in 1:n_cells) {
        
        # Compute differences to every other cell
        di <- sweep(transposed_embedding, 2, transposed_embedding[, i], "-")
        
        # Normalize and scale the differences (L2 norm)
        di <- sweep(di, 2, sqrt(colSums(di^2)), "/") * scale
        di[, i] <- zero_vec  # No distance to itself
        
        # Calculate displacement vector
        transition_shift <- di %*% transition_matrix[, i]
        binarized_shift <- di %*% tpb[, i]

        # Print to debug
        print(paste("Cell:", i))
        print("transition_shift:")
        print(transition_shift)
        print("binarized_shift:")
        print(binarized_shift)
            
        # Final delta embedding
        delta_embedding[i, ] <- transition_shift - binarized_shift
    }
    
    return(delta_embedding)
}


#' Compute Local Coherence of a Vector Field
#'
#' This function calculates the local coherence of a cell-wise vector field
#' (e.g., perturbation vectors) based on a cell–cell neighborhood graph.
#' For each cell, it measures how similar the vectors of its neighbors are
#' to the cell's own vector, using cosine similarity. The result is a
#' per-cell numeric vector of coherence values.
#'
#' @param seurat_obj A \code{Seurat} object containing the neighborhood graph
#'   and the perturbation vectors.
#' @param perturbation_name Name of the perturbation assay or matrix used to
#'   compute the cell-wise vectors (passed to \code{PerturbationVectors}).
#' @param reduction Character. Dimensionality reduction to use for the vector
#'   representation (default: \code{"umap"}).
#' @param graph Character. Name of the graph slot in the Seurat object (e.g.
#'   \code{"RNA_nn"} or \code{"RNA_snn"}) used to define cell neighborhoods.
#' @param n_threads Integer. Number of threads to use (currently not implemented
#'   for parallelization, but reserved for future use).
#' @param arrow_scale Numeric. Scaling factor for arrow (vector) lengths passed
#'   to \code{PerturbationVectors}.
#' @param max_pct Numeric. Maximum percentile of vector magnitudes to keep when
#'   scaling (passed to \code{PerturbationVectors}).
#' @param weighted Logical. If \code{TRUE}, coherence is computed as a
#'   weighted average of neighbor similarities using the edge weights from
#'   the graph; if \code{FALSE}, an unweighted mean is used.
#'
#' @return A numeric vector of length equal to the number of cells, containing
#'   the local coherence score for each cell. Higher values indicate that a
#'   cell’s vector is more aligned with its neighbors’ vectors.
#'
#' @details
#' Local coherence quantifies directional agreement of perturbation vectors
#' within the local neighborhood of each cell. Cosine similarity is used to
#' measure alignment between a cell's vector and those of its neighbors.
#' Weighted coherence accounts for edge weights in the neighborhood graph
#' (e.g. from SNN graphs).
VectorFieldCoherence <- function(
    seurat_obj,
    perturbation_name,
    reduction = 'umap',
    graph = 'RNA_nn',
    n_threads = 4,
    arrow_scale = 1,
    max_pct = 0.90,
    weighted = FALSE
){

    # --- Checks ---
    if (!inherits(seurat_obj, "Seurat")) {
        stop("`seurat_obj` must be a Seurat object.")
    }

    # check reduction
    if (!reduction %in% names(seurat_obj@reductions)) {
        stop(sprintf("Reduction '%s' not found in Seurat object. Available reductions: %s",
                     reduction, paste(names(seurat_obj@reductions), collapse = ", ")))
    }

    # check perturbation_name (assay or data)
    graph_name <- paste0(perturbation_name, "_tp")
    if (!graph_name %in% names(seurat_obj@assays) && 
        !perturbation_name %in% names(seurat_obj@assays)) {
        stop(sprintf("Perturbation '%s' not found in Seurat object assays. Available assays: %s",
                     perturbation_name, paste(names(seurat_obj@assays), collapse = ", ")))
    }

    # check graph slot
    if (!graph %in% names(seurat_obj@graphs)) {
        stop(sprintf("Graph '%s' not found in Seurat object graphs. Available graphs: %s",
                     graph, paste(names(seurat_obj@graphs), collapse = ", ")))
    }

    # get graph
    cell_graph <- Graphs(seurat_obj, slot = graph)

    # neighbor indices and weights
    knn_idx <- lapply(1:nrow(cell_graph), function(x) {
        which(cell_graph[x, ] > 0)
    })
    weights_list <- lapply(1:nrow(cell_graph), function(x) {
        cell_graph[x, cell_graph[x, ] > 0]
    })

    # get perturbation vectors
    vectors <- PerturbationVectors(
        seurat_obj,
        perturbation_name = perturbation_name,
        reduction = reduction, 
        arrow_scale = arrow_scale,
        max_pct = max_pct
    )
    vector_field <- vectors$arsd
    
    # precompute norms
    norms <- sqrt(rowSums(vector_field^2))
    n <- nrow(vector_field)
    coherence <- numeric(n)

    for (i in seq_len(n)) {
        vi <- vector_field[i, , drop = FALSE]
        vi_norm <- norms[i]

        neighbors <- knn_idx[[i]]
        if (length(neighbors) == 0 || vi_norm == 0) {
            coherence[i] <- 0
            next
        }

        vj <- vector_field[neighbors, , drop = FALSE]
        vj_norms <- norms[neighbors]

        # cosine similarities
        dots <- as.matrix(vj) %*% t(as.matrix(vi))  # (k x 1)
        sims <- dots / (vi_norm * vj_norms)

        if (weighted) {
            # use graph weights
            w <- weights_list[[i]]
            w <- w / sum(w)
            coherence[i] <- sum(w * sims, na.rm = TRUE)
        } else {
            coherence[i] <- mean(sims, na.rm = TRUE)
        }
    }

    return(coherence)
}

embArrows_velocyto <- function(emb, tp, arrowScale = 1.0, nthreads = 1L) {
    .Call('_velocyto_R_embArrows', PACKAGE = 'velocyto.R', emb, tp, arrowScale, nthreads)
}



