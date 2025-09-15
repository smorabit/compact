
#' ApplyPerturbation
#'
#' This function applies an in-silico perturbation to selected features in a Seurat object.
#' 
#' @return A dgCMatrix object containing the updated expression matrix with the applied perturbations
#' 
#' @param seurat_obj A Seurat object
#' @param exp A features by cells matrix containing the observed expression matrix.
#' @param features The selected features to apply the perturbation on.
#' @param perturb_dir A numeric determining the type of perturbation to apply. Negative values for knock-down, positive for knock-in, and 0 for knock-out.
#' @param slot Slot to extract data for aggregation. Default = 'counts'
#' @param assay Assay in seurat_obj containing expression information.
#' 
#' @details 
#'
#' @import Seurat
#' @export
ApplyPerturbation <- function(
    seurat_obj,
    exp,
    features,
    perturb_dir,
    cells_use=NULL,
    group.by = NULL,
   # group_name = NULL,
    layer = 'counts',
    slot = 'counts',
    assay = 'RNA'
){

   # if we are doing a knock-out:
    if(perturb_dir == 0){
        exp[features,] <- 0
        return(exp)
    }

    # which cells are we using?
    # if(is.null(group.by)){
    #   cells_use <- colnames(seurat_obj)
    # } else{
    #     cells_use <- seurat_obj@meta.data %>% 
    #       subset(get(group.by) %in% group_name) %>%
    #       rownames
    # }

    # split the exp matrix based on selected features
    exp_hubs <- exp[features,cells_use]
    exp_non <- exp[!(rownames(exp) %in% features),cells_use]

    # # define groups based on group.by
    # if(is.null(group.by)){
    #     group.by <- 'fake_group'
    #     seurat_obj@meta.data[,group.by] <- "all"
    #     groups <- c("all")
    # } else{
         groups <- unique(as.character(seurat_obj@meta.data[,group.by]))
    # }

    # initialize progress bar
    pb <- utils::txtProgressBar(min = 0, max = length(features), style = 3, width = 50, char = "=")

    # for each feature, model expression 
    # TODO: add options for other models aside from ZINB like Poisson etc.
    # TODO: add option to return the model (just in case?)
    # TODO: should the SampleZINB function return a sparse vector?
    sim_data <- lapply(1:length(features), function(i){
        
        feature <- features[i]
        print(feature)

        # initialize 
        ysim <- c()
        cell_names <- c()

        # get the current group:
        for(cur_group in groups){

            print(cur_group)

            cur_cells <- subset(seurat_obj@meta.data, get(group.by) == cur_group) %>% rownames()

            # get the obserbed gene expression
            if(length(features) == 1){
                yobs <- exp_hubs[cur_cells]
            } else{
                yobs <- exp_hubs[feature,cur_cells]
            }

            # model expression as a ZINB
            model <- ModelZINB(
                seurat_obj, 
                feature=feature, cells_use=cur_cells, 
                layer=layer, slot=slot
            )
            
            # simulate data by samping the distribution 
            cur_ysim <- SampleZINB(model, yobs)
            cur_ysim <- cur_ysim[1:length(cur_cells)]
            
            # this line doesn't work
            ysim <- c(ysim, cur_ysim)
            cell_names <- c(cell_names, cur_cells)
            print(length(ysim))
        }

        # set the names for ysim 
        names(ysim) <- cell_names
        
        # update progress bar
        setTxtProgressBar(pb, i)

        # return the simulated data
        ysim
    })

    # close progress bar
    close(pb)

    # is there one gene or more than 1 being perturbed?
    if(length(sim_data) > 1){
        sim_data <- do.call(rbind, sim_data)

        # convert to sparse matrix:
        sim_data <- Matrix::Matrix(sim_data, sparse=TRUE)

        # re-order cells to match the original matrix:
        sim_data <- sim_data[,colnames(exp_hubs)]

        # is this a knock-in or knock-down?
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
