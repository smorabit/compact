

#' ModelZINB
#'
#' This function models the expression of a selected feature as a zero-inflated negative binomial (ZINB) distribution.
#' 
#' @return The zero-inflated negative binomial model fit to the observed data from the selected feature.
#' 
#' @param seurat_obj A Seurat object
#' @param feature The selected feature to model, must be present in rownames(seurat_obj)
#' @param feature List of cell barcodes from the Seurat object
#' @param slot Slot to extract data for aggregation. Default = 'counts'
#' 
#' @details 
#' ModelZINB is an internal helper function that calls on zeroinfl from the pscl package to 
#' model expression of a selected feature.
#'
#' @import Seurat
#' @import pscl
#' @export
ModelZINB <- function(
    seurat_obj,
    feature,
    cells_use=NULL,
    layer = 'counts',
    slot = 'counts',
    add_zero = TRUE
){

    # check slot
    if(!(slot %in% c('counts', 'data', 'scale.data'))){
        stop(paste0("Invalid slot (", slot, "). Valid options for slot: counts, data, scale.data "))
    }

    # check that the feature is in the rownames Seurat object:
    if(!(feature %in% rownames(seurat_obj))){
        stop(paste0("Invalid feature, ", feature, " not found in rownames(seurat_obj)."))
    }

    if(is.null(cells_use)){
        cells_use <- colnames(seurat_obj)
    }

    # validate cells_use barcodes
    invalid_cells <- setdiff(cells_use, colnames(seurat_obj))
    if(length(invalid_cells) > 0){
        preview <- paste(head(invalid_cells, 5), collapse = ', ')
        suffix  <- if(length(invalid_cells) > 5) paste0(" ... (", length(invalid_cells) - 5, " more)") else ""
        stop(paste0("cells_use contains barcodes not found in colnames(seurat_obj): ", preview, suffix))
    }

    # get original expression matrix
    if(hdWGCNA::CheckSeurat5()){
        cur_x <- FetchData(seurat_obj, feature, layer=layer, cells=cells_use)[,1]
    } else{
        cur_x <- FetchData(seurat_obj, feature, slot=slot, cells=cells_use)[,1]
    }

    if(add_zero){
        cur_x <- c(cur_x, 0)
    }

    # ZINB cannot be fit to all-zero data; catch this before zeroinfl throws
    # an opaque convergence error
    if(sum(cur_x) == 0){
        stop(paste0("Cannot fit ZINB model: all expression values for '", feature,
                    "' are zero in the selected cells. ",
                    "This gene may not be expressed in this cell group."))
    }

    # fit data to zero-inflated negative binomial distribution
    zinb_model <- tryCatch(
        pscl::zeroinfl(x ~ . | 1, dist = 'negbin', data = data.frame(x = cur_x)),
        error = function(e){
            stop(paste0("ZINB model fitting failed for gene '", feature, "': ",
                        conditionMessage(e)))
        }
    )

    # return the model:
    zinb_model

}

#' SampleZINB
#'
#' This function simulated expression values based on a ZINB model that has been fit to gene expression data.
#' 
#' @return A numeric vector of simulated expression data based on a ZINB model
#' 
#' @param model a model object
#' @param ncells number of cells to simulate expression values for. By default will simulate data for each cell in the input dataset.
#'
#' @import VGAM
#' @import pscl
#' @export
SampleZINB <- function(
    model,
    yobs,
    ncells = NULL
){

    if(is.null(ncells)){
        ncells <- model$n
    }

    # rzinegbin requires munb > 0; return zeros directly for unexpressed genes
    if(mean(yobs) == 0){
        return(rep(0L, ncells))
    }

    # get the parameters for simulating data from this model
    theta <- model$theta
    zero_intercept <- plogis(model$coefficients$zero)

    # simulating data based on the zinb model fit
    ysim <- VGAM::rzinegbin(
        n = ncells,
        munb = mean(yobs),
        size = theta,
        pstr0 = zero_intercept
    )

    # return the simulated data
    ysim

}
