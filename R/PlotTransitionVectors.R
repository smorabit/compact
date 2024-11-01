
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
    point_size = 1

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
    emb <- Reductions(seurat_obj, reduction)@cell.embeddings
    colnames(emb) <- paste0('emb_', 1:ncol(emb))

    # make a grid of arrows 
    grid.df <- grid_vectors(emb[rownames(ars),], arsd, resolution=grid_resolution)

    # make a dataframe to plot with ggplot
    plot_df <- as.data.frame(emb)

    # make a dataframe to plot with ggplot
    plot_df <- as.data.frame(emb)
    plot_df$value <- seurat_obj@meta.data[,color.by]

    # plot the scatter plot
    p <- ggplot(plot_df, aes(x=emb_1, y=emb_2, color=value)) +
        geom_point(alpha=point_alpha, size=point_size) 
    
    # add the arrows
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
        size=0.25
    )
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
    emb <- Reductions(seurat_obj, reduction)@cell.embeddings

    # get the graph from the seurat object
    tp <- Graphs(seurat_obj, graph_name)

    # run the helper function
    arsd <- data.frame(t(embArrows_velocyto(
        emb,
        tp,
        arrow_scale,
        n_threads
    )))

    # format the dataframes
    rownames(arsd) <- rownames(emb)
    ars <- data.frame(cbind(emb,emb+arsd))
    colnames(ars) <- c('x0','y0','x1','y1')
    colnames(arsd) <- c('xd','yd')
    rownames(ars) <- rownames(emb)


    max_vect <- quantile(abs(as.matrix(arsd)), max_pct)
    arsd[arsd > max_vect] <- max_vect
    arsd[arsd < -max_vect] <- -max_vect
    arsd <- arsd / max_vect  

    # exclude NAs:
    ars <- na.omit(ars)
    arsd <- na.omit(arsd)

    # return 
    list(ars=ars, arsd=arsd)

}




#' @export
#' @importFrom S4Vectors selfmatch DataFrame
#' @importFrom stats median
grid_vectors <- function(x, embedded, resolution=40, scale=TRUE, as.data.frame=TRUE, return.intermediates=FALSE) {
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

    if (return.intermediates) {
        if (as.data.frame) {
            warning("ignoring 'as.data.frame=TRUE' for 'return.intermediates=TRUE'")
        }
        return(list(start=pos, end=pos + vec, limits=limits, delta=delta,
                    categories=categories, grp=grp, vec=vec))
    } else {
        FUN <- if (as.data.frame) data.frame else list
        return(FUN(start=pos, end=pos + vec))
    }
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

embArrows_velocyto <- function(emb, tp, arrowScale = 1.0, nthreads = 1L) {
    .Call('_velocyto_R_embArrows', PACKAGE = 'velocyto.R', emb, tp, arrowScale, nthreads)
}