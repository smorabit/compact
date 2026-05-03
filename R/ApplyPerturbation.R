
#' Apply In-Silico Perturbation
#'
#' This function applies an in-silico perturbation (knock-out, knock-down, or knock-in) 
#' to selected features in a Seurat object.
#' 
#' @param seurat_obj A Seurat object containing the dataset.
#' @param exp A features-by-cells matrix (typically a sparse matrix) containing the observed expression data.
#' @param features Character vector. The selected features to apply the perturbation on.
#' @param perturb_dir Numeric. Determines the type of perturbation to apply. Negative values for knock-down, positive for knock-in, and 0 for knock-out.
#' @param cells_use Character vector. Specific cells to apply the perturbation to. If \code{NULL}, defaults to handling internally.
#' @param group.by Character. Column in \code{seurat_obj@meta.data} used to group cells for modeling.
#' @param layer Character. The Seurat v5 layer to extract data from. Default is \code{'counts'}.
#' @param slot Character. The Seurat v4 slot to extract data from. Default is \code{'counts'}.
#' @param assay Character. The assay in \code{seurat_obj} containing expression information. Default is \code{'RNA'}.
#' @param n_workers Integer. Number of parallel workers for fitting and sampling the ZINB model across
#' features. Uses fork-based parallelism (\code{parallel::mclapply}) so memory usage does not scale with
#' worker count. Default is \code{1} (serial). Values > 1 disable the progress bar.
#'
#' @details
#' The function models the baseline expression of target features using a Zero-Inflated
#' Negative Binomial (ZINB) distribution to capture dropout and overdispersion characteristics
#' typical of single-cell RNA-seq data .
#' It samples from this distribution, scales the values by \code{perturb_dir}, and updates the
#' original expression matrix. Any perturbed counts that fall below zero are strictly bounded to zero.
#'
#' @return A \code{dgCMatrix} object containing the updated expression matrix with the applied perturbations.
#'
#' @import Seurat
#' @importFrom Matrix Matrix
#' @importFrom parallel mclapply
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
ApplyPerturbation <- function(
    seurat_obj,
    exp,
    features,
    perturb_dir,
    cells_use=NULL,
    group.by = NULL,
    layer = 'counts',
    slot = 'counts',
    assay = 'RNA',
    n_workers = 1
){

   # if we are doing a knock-out:
    if(perturb_dir == 0){
        exp[features,] <- 0
        return(exp)
    }

    # default to all cells if not specified
    if(is.null(cells_use)) cells_use <- colnames(seurat_obj)

    # default group.by to a single group covering all cells
    if(is.null(group.by)){
        group.by <- 'fake_group'
        seurat_obj@meta.data[, group.by] <- 'all'
    }

    # split the exp matrix based on selected features
    exp_hubs <- exp[features,cells_use]
    exp_non <- exp[!(rownames(exp) %in% features),cells_use]

    groups <- unique(as.character(seurat_obj@meta.data[,group.by]))

    # per-feature ZINB fit + sample; each feature is independent so this is
    # embarrassingly parallel. fork-based mclapply (n_workers > 1) keeps memory
    # flat because workers share the parent's address space via copy-on-write
    # and only read seurat_obj / exp_hubs — no large objects are copied.
    # TODO: add options for other models aside from ZINB like Poisson etc.
    # TODO: add option to return the model (just in case?)
    # TODO: should the SampleZINB function return a sparse vector?
    worker <- function(i){
        feature <- features[i]

        ysim       <- c()
        cell_names <- c()

        for(cur_group in groups){
            cur_cells <- subset(seurat_obj@meta.data, get(group.by) == cur_group) %>% rownames()

            yobs <- if(length(features) == 1) exp_hubs[cur_cells] else exp_hubs[feature, cur_cells]

            model <- ModelZINB(
                seurat_obj,
                feature   = feature,
                cells_use = cur_cells,
                layer     = layer,
                slot      = slot
            )

            # ncells is passed explicitly because model$n includes the
            # artificial zero appended by ModelZINB (add_zero = TRUE)
            cur_ysim   <- SampleZINB(model, yobs, ncells = length(cur_cells))
            ysim       <- c(ysim, cur_ysim)
            cell_names <- c(cell_names, cur_cells)
        }

        names(ysim) <- cell_names
        ysim
    }

    if(n_workers == 1){
        pb       <- utils::txtProgressBar(min = 0, max = length(features), style = 3, width = 50, char = "=")
        sim_data <- lapply(seq_along(features), function(i){
            result <- worker(i)
            utils::setTxtProgressBar(pb, i)
            result
        })
        close(pb)
    } else {
        sim_data <- parallel::mclapply(
            seq_along(features), worker,
            mc.cores      = n_workers,
            mc.preschedule = FALSE   # avoids load imbalance since ZINB fit time varies by gene
        )
    }

    # is there one gene or more than 1 being perturbed?
    if(length(sim_data) > 1){
        sim_data <- do.call(rbind, sim_data)

        # convert to sparse matrix:
        sim_data <- Matrix::Matrix(sim_data, sparse=TRUE)

        # re-order cells to match the original matrix:
        sim_data <- sim_data[,colnames(exp_hubs)]

        # is this a knock-in or knock-down?
        #
        # <TODO> we need to handle the case where perturb_dir is a decimal
        # just need to round the results so it remains a valid counts matrix
        sim_data <- sim_data * perturb_dir

        # apply the perturbation to the expression matrix
        delta_hub <- exp_hubs + sim_data

        # if in a knock-down experiment any feature fell below 0 expression, make it 0.
        delta_hub[delta_hub < 0 ] <- 0

        # update the expression matrix with the perturbed values
        exp_new <- rbind(delta_hub, exp_non)

    } else{

        sim_data <- sim_data[[1]]

        # re-order cells to match the original matrix
        sim_data <- sim_data[names(exp_hubs)]

        # is this a knock-in or knock-down?
        sim_data <- sim_data * perturb_dir

        # apply the perturbation to the expression matrix
        delta_hub <- exp_hubs + sim_data

        # if in a knock-down experiment any feature fell below 0 expression, make it 0.
        delta_hub[delta_hub < 0 ] <- 0

        # update the expression matrix with the perturbed values
        exp_new <- rbind(delta_hub, exp_non)
        rownames(exp_new)[1] <- features
    }

    #get the data for the cells that we did not perturb
    if(!all(colnames(seurat_obj) %in% cells_use)){
      exp_other <- exp[,!(colnames(exp) %in% cells_use)]
      exp_new <- cbind(exp_new, exp_other)
    }

    # re-order the genes to match the original matrix
    exp_new <- exp_new[rownames(exp),colnames(seurat_obj)]
    exp_new <- exp_new[rownames(exp),]

    # return the expression matrix where the perturbation has been applied: 
    exp_new
}
