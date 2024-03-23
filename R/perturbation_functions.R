
#' ModulePerturbation
#'
#' This function enables in-silico gene expression perturbation analysis using 
#' a co-expression network. 
#' 
#' @return A Seurat object containing the in-silico perturbation results as a new assay.
#' 
#' @param seurat_obj A Seurat object
#' @param mod Name of the co-expression module to perturb
#' @param perturb_dir A numeric determining the type of perturbation to apply. Negative values for knock-down, positive for knock-in, and 0 for knock-out.
#' @param perturbation_name A name for the in-silico perturbation that will be stored in the Seurat obejct
#' @param graph Name of the cell-cell graph in the Graphs(seurat_obj)
#' @param group.by A string containing the name of a column in the Seurat object with cell groups (clusters, cell types, etc).
#' @param group_name A string containing a group present in the provided group.by column. A character vector can be provided to select multiple groups at a time.
#' @param n_hubs The number of hub genes to perturb from the selected co-expression module.
#' @param n_iters The number of times to apply the signal propagation.
#' @param corr_sigma A numeric scaling factor for the correlation matrix.
#' @param n_threads Number of threads for the correlation calculation
#' @param slot Slot to extract data for aggregation. Default = 'counts'
#' @param assay Assay in seurat_obj containing expression information.
#' @param return_delta Logical indicating whether or not to return the "Delta" matrix in the Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' 
#' @details 
#' Following co-expression network analysis with hdWGCNA, ModulePerturbation allows us to perform 
#' in-silico gene expression perturbation analysis. This analysis consists of several key steps.
#' 
#' 1. Apply a primary perturbation to the hub genes of a selected module. In this step, 
#' we model the observed gene expression of the hub genes using a zero-inflated negative binomial (ZINB),
#' or other distributions. We then simulate new exptession data by sampling this distribution.
#' The simulated expression matrix is then multiplied by the perturb_dir, and the final perturbation expression 
#' matrix is computed by adding the observed and simulated expression matrices.
#' 
#' 2. Apply a secondary perturbation to the rest of the co-expression module using a signal propagation
#' algorithm. Given the gene-gene co-expression network from hdWGCNA, we can propagate the perturbation signal 
#' throughout the network by computing the dot product between the network adjacency matrix and the perturbation 
#' expression matrix. This can be performed over several iterations using n_iters.
#' 
#' 3. Compute the cell-cell transition probabilities. 
#'
#' @import Seurat
#' @export 
#' ModulePerturbation
ModulePerturbation <- function(
    seurat_obj,
    mod,
    perturb_dir,
    perturbation_name,
    graph,
    group.by = NULL,
    group_name = NULL,
    n_hubs = 5,
    n_iters = 3,
    delta_scale = 0.2,
    corr_sigma=0.05,
    n_threads=4,
    slot = 'counts',
    assay = 'RNA',
    #all_features=FALSE,
    return_delta = FALSE, # TODO: this does nothing for now
    wgcna_name=NULL
){

    # set as active assay if wgcna_name is not given
    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

    # check assay
    if(!(assay %in% names(seurat_obj@assays))){
        stop(paste0("Invalid assay (", assay, "). Assays in Seurat object: ", paste(names(seurat_obj@assays), collapse=', ')))
    }

    # switch to this assay:
    orig_assay <- DefaultAssay(seurat_obj)
    DefaultAssay(seurat_obj) <- assay

    # check slot
    if(!(slot %in% c('counts', 'data', 'scale.data'))){
        stop(paste0("Invalid slot (", slot, "). Valid options for slot: counts, data, scale.data "))
    }

    # check perturb_dir 
    if(!is.numeric(perturb_dir)){
        stop(paste0('Invalid choice for perturb_dir. Valid choices are positive numbers (knock-in), negative numbers (knock-down), or 0 (knock-out).'))
    }


    # checks 
    # TODO: which checks should go here and which should go inside the functions?
    # it would be nice to do them all here so it errors out immediately.

    # get the TOM:
    TOM <- GetTOM(seurat_obj, wgcna_name)

    # get modules 
    modules <- GetModules(seurat_obj, wgcna_name)

    # get top n hub genes in selected module 
    hubs <- GetHubGenes(seurat_obj, n=n_hubs, wgcna_name=wgcna_name) %>% 
        subset(module == mod)
    hub_genes <- hubs$gene_name

    # get the non-hubs in this module:
    non_hub_genes <- subset(modules, module == mod & !(gene_name %in% hub_genes)) %>% .$gene_name
    module_genes <- subset(modules, module == mod) %>% .$gene_name

    # which cells are we selecting to apply the perturbation?
    if(is.null(group.by)){
      cells_use <- colnames(seurat_obj)
    } else{
        cells_use <- seurat_obj@meta.data %>% 
          subset(get(group.by) %in% group_name) %>%
          rownames
    }

    ###########################################################################
    # Set up the observed expression matrix
    # TODO: make this into its own function
    ###########################################################################

    # get the expression matrix:
    exp <- GetAssayData(seurat_obj, slot=slot, assay = assay)

    #print(mean(exp['PTPN2',]))


    ###########################################################################
    # Part 1: apply the perturbation to the selected hub genes:
    ###########################################################################

    print('Applying in-silico perturbation to hub genes...')
    exp_per <- ApplyPerturbation(
        seurat_obj,
        exp,
        features = hub_genes,
        perturb_dir = perturb_dir,
        cells_use = cells_use,
       # group.by = group.by,
       # group_name = group_name,
        slot = slot,
        assay = assay
    )

   # print(mean(exp_per['PTPN2',]))


    ###########################################################################
    # Part 2: apply signal propagation throughout this module 
    ###########################################################################

    # get the TOM for the genes in this module
    cur_TOM <- TOM[module_genes, module_genes]

    # set the TOM to zero for non hub genes so only the connections with the
    # perturbed hub genes play a role
    #cur_TOM[,non_hub_genes] <- 0
   # cur_TOM[hub_genes,] <- 0

    print('Applying signal propagation throughout co-expression network...')
    exp_prop <- ApplyPropagation(
        seurat_obj,
        exp[module_genes,cells_use],
        exp_per[module_genes,cells_use],
        network = cur_TOM,
        perturb_dir = perturb_dir,
        delta_scale = delta_scale,
        n_iters = n_iters
    )

    exp_prop2 <- exp_prop
  #  print(mean(exp_prop['PTPN2',]))

    # something messed up happens after here

    if(!all(colnames(seurat_obj) %in% cells_use)){
        exp_prop_other <- exp[module_genes,setdiff(colnames(seurat_obj), cells_use)]
        exp_prop <- cbind(exp_prop, exp_prop_other)
        exp_prop <- exp_prop[,colnames(seurat_obj)]
    }

    # append the expression matrices:
    exp_simulated <- rbind(
        exp_per[!(rownames(exp_per) %in% module_genes),], # genes that aren't in this module
        exp_prop # genes from this module with perturbations
    )

    # make sure the order matches the original expression matrix
    exp_simulated <- exp_simulated[rownames(seurat_obj),colnames(seurat_obj)]

    # add perturbation assay to the Seurat object:
    perturb_assay <- CreateAssayObject(
        exp_simulated,
        assay = perturbation_name
    )
    seurat_obj[[perturbation_name]] <- perturb_assay

    # normalize the perturbation assay
    seurat_obj <- NormalizeData(seurat_obj, perturbation_name)

    ###########################################################################
    # Part 3: compute transition probabilities
    ###########################################################################

    # print('Computing cell-cell transition probabilities based on the perturbation...')
    seurat_obj <- PerturbationTransitions(
        seurat_obj,
        perturbation_name,
        features=module_genes,
        graph=graph, 
        #all_features=all_features,
        corr_sigma=corr_sigma,
        n_threads=n_threads,
        slot='data', # should we make this an option?
        assay=assay
    )

    # return the Seurat object
    seurat_obj

}


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
ApplyPerturbation <- function(
    seurat_obj,
    exp,
    features,
    perturb_dir,
    cells_use=NULL,
   # group.by = NULL,
   # group_name = NULL,
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

    # initialize progress bar
    pb <- utils::txtProgressBar(min = 0, max = length(features), style = 3, width = 50, char = "=")

    # for each feature, model expression 
    # TODO: add options for other models aside from ZINB like Poisson etc.
    # TODO: add option to return the model (just in case?)
    # TODO: should the SampleZINB function return a sparse vector?
    sim_data <- lapply(1:length(features), function(i){
        
        feature <- features[i]

        # get the obserbed gene expression
        if(length(features) == 1){
            yobs <- exp_hubs
        } else{
            yobs <- exp_hubs[feature,]
        }

        # model expression as a ZINB
        model <- ModelZINB(seurat_obj, feature=feature, cells_use=cells_use, slot=slot)
        
        # simulate data by samping the distribution 
        ysim <- SampleZINB(model, yobs)

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


#' ApplyPropagation
#'
#' This function applies an in-silico perturbation to selected features in a Seurat object.
#' 
#' @return A dgCMatrix object containing the updated expression matrix with the applied perturbations
#' 
#' @param seurat_obj A Seurat object
#' @param exp A features by cells matrix containing the observed expression matrix.
#' @param exp_per A features by cells matrix containing the perturbation results from ApplyPerturbation.
#' @param network A gene-gene network to apply the signal propagation.
#' @param n_iters The number of times to apply the signal propagation.
#' 
#' @details 
#'
#' @import Seurat
ApplyPropagation <- function(
    seurat_obj,
    exp,
    exp_per,
    network,
    perturb_dir = perturb_dir,
    n_iters = 3,
    delta_scale = 0.2
){

    # todo: some checks for the network?
    # check that the network has the right genes

    # compute the difference between the perturbation and the observed expression
    delta <- exp_per - exp

    # backups
    delta_orig <- delta 
    exp_orig <- exp 
    exp_per_orig <- exp_per

    # what is the sign of the perturbation dir?
    perturb_sign <- sign(perturb_dir)

    # run the signaling processing step iteratively
    for(i in 1:n_iters){  

        # compute the dot product between the TOM coefficients and the exp matrix
        delta <- network %*% delta

        # penalize the delta, or else the values rapidly get too large
        delta <- delta * delta_scale

        # update the expression matrix:
        exp_update <- exp_per + delta # (delta * perturb_sign)
        exp_update[exp_update < 0] <- 0

        exp_update <- round(exp_update)

        # update the delta
        delta <- exp_update - exp_per

    }
    
    # return the matrix with the signal propagation applied
    exp_prop <- exp_per + delta
    exp_prop[exp_prop < 0] <- 0
    exp_prop
}

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
ModelZINB <- function(
    seurat_obj,
    feature,
    cells_use=NULL,
    slot = 'counts'
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

    # get the expression profile for this feature
    cur_x <- FetchData(seurat_obj, feature , slot=slot, cells=cells_use)[,1]
    
    # fit data to zero-inflated negative binomial distribution
    zinb_model <- pscl::zeroinfl(cur_x ~ rep(0, length(cur_x)) | 1, dist='negbin')

    # return the model:
    zinb_model

}

#' SampleZINB
#'
#' This function simulated expression values based on a ZINB model that has been fit to gene expression data.
#' 
#' @return A numeric vector of simulated expression data based on a ZINB model
#' 
#' @param model
#' @param ncells number of cells to simulate expression values for. By default will simulate data for each cell in the input dataset.
#'
#' @import VGAM
#' @import pscl
SampleZINB <- function(
    model,
    yobs,
    ncells = NULL
){

    if(is.null(ncells)){
        ncells <- model$n
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



#' PerturbationTransitions
#'
#' This function computes cell-cell transition probabilities based on an in-silico perturbation experiment.
#' 
#' @return A Seurat object with the cell-cell transition probabilities stored in the Graphs slot of the Seurat object.
#' 
#' @param seurat_obj A Seurat object
#' @param perturbation_name A name for the in-silico perturbation that will be stored in the Seurat obejct
#' @param features Selected features to use for the transition probability calculation
#' @param graph Name of the cell-cell graph in the Graphs(seurat_obj)
#' @param corr_sigma A numeric scaling factor for the correlation matrix.
#' @param n_threads Number of threads for the correlation calculation
#' @param slot Slot to extract data for aggregation. Default = 'data'
#' @param assay Assay in seurat_obj containing expression information.
#'
#' @import Matrix
#' @import Seurat
PerturbationTransitions <- function(
    seurat_obj,
    perturbation_name,
    features,
    graph, # graph in the seurat object
   # all_features = FALSE,
    corr_sigma=0.05,
    n_threads=4,
    slot='data',
    assay="RNA"
){

    # check for selected Graph 
    cell_graph <- Graphs(seurat_obj, slot=graph)
    cell_graph <- Matrix::Matrix(cell_graph)
    diag(cell_graph) <- 1

    # if(all_features){
    #     features <- rownames(seurat_obj)
    # }
   
    # get the observed and the perturbed expression matrices:
    exp_obs <- GetAssayData(seurat_obj, assay=assay, slot=slot)
    exp_per <- GetAssayData(seurat_obj, assay=perturbation_name, slot=slot)

    # subset by selected features
    # and convert to dense matrix (TODO: figure out a way to avoid this?)
    exp_obs <- as.matrix(exp_obs[features,])
    exp_per <- as.matrix(exp_per[features,])

    delta <- exp_per - exp_obs

    # run the Velocyto colDeltaCor function
    cc <- colDeltaCor(exp_obs, delta, nthreads=n_threads)

    # fill the diagnoal with zeros (cells won't transition to self)
    diag(cc) <- 0

    # compute transition probs between cells
    tp <- exp(cc / corr_sigma) * cell_graph
    tp <- t(t(tp)/Matrix::colSums(tp)); # tp shows transition from a given column cell to different row cells
    tp <- as(tp,'dgCMatrix') # cast to sparse matrix

    # add rownames / colnames 
    rownames(tp) <- colnames(seurat_obj)
    colnames(tp) <- colnames(seurat_obj)

    # convert to Seurat Graph object, and add it to the Seurat object with the perturbation name
    tp <- Seurat::as.Graph(tp)
    graph_name <- paste0(perturbation_name, '_tp')
    seurat_obj@graphs[graph_name] <- tp

    # return the Seurat object 
    seurat_obj

}



PerturbationVectors <- function(
    seurat_obj,
    perturbation_name,
    reduction = 'umap',
    n_threads=4,
    arrow_scale = 1
){

    # TODO: check reduction 
    graph_name <- paste0(perturbation_name, '_tp')

    # get the 2D embedding
    emb <- Reductions(seurat_obj, reduction)@cell.embeddings

    # get the graph from the seurat object
    tp <- Graphs(seurat_obj, graph_name)

    # run the helper function
    arsd <- data.frame(t(embArrows(
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

    # exclude NAs:
    ars <- na.omit(ars)
    arsd <- na.omit(arsd)

    # return 
    list(ars=ars, arsd=arsd)

}


PlotPerturbationVectors <- function(
    seurat_obj,
    perturbation_name,
    reduction = 'umap',
    n_threads=4,
    arrow_scale = 1,
    grid_n = 25,
    grid_sd = NULL,
    min_arrow_size = NULL,
    min_grid_arrow_length = NULL,
    min_grid_cell_mass = 1,
    arrow_lwd = 1
){

    # maybe move this outside of this function so we don't have to run it every time we want to plot 
    vectors <- PerturbationVectors(
        seurat_obj,
        perturbation_name = 'perturb'
    )
    ars <- vectors$ars 
    arsd <- vectors$arsd

    ################################################################################
    # create arrow grid
    ################################################################################

    rx <- range(c(range(ars$x0),range(ars$x1)))
    ry <- range(c(range(ars$y0),range(ars$y1)))
    gx <- seq(rx[1],rx[2],length.out=grid.n)
    gy <- seq(ry[1],ry[2],length.out=grid.n)

    # for each grid point calculate Gaussian-weighted delta average
    if(is.null(grid.sd)) {
    grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
    #  cat("grid.sd=",grid.sd," ")
    }

    if(is.null(min.arrow.size)) {
    min.arrow.size <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)*1e-2;
    cat("min.arrow.size=",min.arrow.size," ")
    }

    if(is.null(max.grid.arrow.length)) {
    max.grid.arrow.length <- sqrt(sum((par('pin')/c(length(gx),length(gy)))^2))*0.25
    cat("max.grid.arrow.length=",max.grid.arrow.length," ")
    }


    garrows <- do.call(rbind,lapply(gx,function(x) {

    # cell distances (rows:cells, columns: grid points)
    cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
    cw <- dnorm(cd,sd=grid.sd)

    # calculate x and y delta expectations
    gw <- Matrix::colSums(cw)
    cws <- pmax(1,Matrix::colSums(cw));
    gxd <- Matrix::colSums(cw*arsd$xd)/cws
    gyd <- Matrix::colSums(cw*arsd$yd)/cws

    al <- sqrt(gxd^2+gyd^2);
    vg <- gw>=min.grid.cell.mass & al>=min.arrow.size


    # run this to get all grid 
    #vg <- vg == vg

    # this subsets by valid points
    cbind(rep(x,sum(vg)),gy[vg],x+gxd[vg],gy[vg]+gyd[vg])

    }))
    colnames(garrows) <- c('x0','y0','x1','y1')


    ################################################################################
    # ggplot
    ################################################################################

    plot_df <- as.data.frame(emb)
    plot_df$pseudotime <- seurat_mg$pseudotime

    plot_df <- plot_df[rownames(ars),]

    p <- ggplot(plot_df, aes(x=UMAP_1, y=UMAP_2, color=pseudotime)) +
    geom_point() +
    scale_color_gradientn(colors=plasma(256))
    g <- ggplot_build(p)
    g_df <- g$data[[1]]

    plot_df$color <- g_df$colour
    cell.colors <- plot_df$color; names(cell.colors) <- rownames(plot_df)




    # compute arrow lengths
    alen <- pmin(max.grid.arrow.length,sqrt( ((garrows[,3]-garrows[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + ((garrows[,4]-garrows[,2])*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2))


    arrow_df <- as.data.frame(garrows)
    arrow_df$length <- alen

    p <- ggplot(plot_df, aes(x=UMAP_1, y=UMAP_2, color=pseudotime)) +
    geom_point(alpha=0.25, size=0.5) +
    scale_color_gradientn(colors=plasma(256)) +
    geom_segment(
        data = arrow_df,
        inherit.aes=FALSE,
        aes(x=x0, y=y0, xend=x1, yend=y1), arrow=grid::arrow(length=unit(0.1, "cm")), size=0.25
    ) +
    coord_equal() +
    umap_theme() +
    NoLegend() +
    ggtitle('MG-M4 knock-down')



}










################################################################################
# Velocyo C++ helper functions
# Is there a way that we could remove the dependency on Velocyto?
################################################################################

colDeltaCor <- function(e, d, nthreads = 1L) {
    .Call('_velocyto_R_colDeltaCor', PACKAGE = 'velocyto.R', e, d, nthreads)
}

embArrows <- function(emb, tp, arrowScale = 1.0, nthreads = 1L) {
    .Call('_velocyto_R_embArrows', PACKAGE = 'velocyto.R', emb, tp, arrowScale, nthreads)
}