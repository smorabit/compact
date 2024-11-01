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
    use_velocyto=TRUE,
    use_graph_tp = FALSE,
    layer = 'counts',
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
    # TODO: make this into its own function (?)
    ###########################################################################

    # get the expression matrix:
    if(hdWGCNA::CheckSeurat5()){
        exp <- Seurat::GetAssayData(seurat_obj, assay = assay, layer=layer)
    } else{
        exp <- Seurat::GetAssayData(seurat_obj, assay = assay, slot=slot)
    }


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
        layer = layer,
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

    print('Computing cell-cell transition probabilities based on the perturbation...')
    seurat_obj <- PerturbationTransitions(
        seurat_obj,
        perturbation_name,
        features=module_genes,
        graph=graph, 
        #all_features=all_features,
        use_velocyto=use_velocyto,
        use_graph_tp = use_graph_tp,
        corr_sigma=corr_sigma,
        n_threads=n_threads,
        layer='data',
        slot='data', 
        assay=assay
    )

    # return the Seurat object
    seurat_obj

}
