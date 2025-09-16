#' TFPerturbation
#'
#' This function enables in-silico transcription factor (TF) perturbation analysis 
#' using a regulatory network derived from co-expression and regulon information.
#' 
#' @return A Seurat object containing the in-silico TF perturbation results as a new assay.
#' 
#' @param seurat_obj A Seurat object.
#' @param selected_tf The name of the transcription factor (TF) to perturb.
#' @param perturb_dir A numeric determining the type of perturbation to apply. Negative values for knock-down, positive for knock-in, and 0 for knock-out.
#' @param perturbation_name A name for the in-silico perturbation that will be stored in the Seurat object.
#' @param graph Name of the cell-cell graph in the Graphs(seurat_obj). Default = "RNA_nn".
#' @param n_iters The number of times to apply the signal propagation throughout the TF regulatory network. Default = 1.
#' @param delta_scale A numeric factor scaling the perturbation during propagation. Default = 1.
#' @param corr_sigma A numeric scaling factor for the correlation matrix. Default = 0.05.
#' @param n_threads Number of threads for the correlation calculation. Default = 4.
#' @param use_velocyto Logical indicating whether to compute velocity-based transition probabilities. Default = TRUE.
#' @param use_graph_tp Logical indicating whether to use the graph topology for transition probabilities. Default = FALSE.
#' @param depth The depth of the regulatory network to use for target identification. Default = 2.
#' @param target_type A string specifying the type of targets to include ("upstream", "downstream", or "both"). Default = "both".
#' @param use_regulons Logical indicating whether to use regulons for TF-target relationships. Default = TRUE.
#' @param layer Layer of the assay containing expression data. Default = "counts".
#' @param slot Slot of the assay containing expression data. Default = "counts".
#' @param assay Assay in seurat_obj containing expression information. Default = "RNA".
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. If NULL, uses the active WGCNA experiment.
#' 
#' @details 
#' TFPerturbation enables in-silico transcription factor perturbation by simulating changes 
#' in TF activity and propagating these changes throughout the associated regulatory network. 
#' The workflow consists of the following steps:
#' 
#' 1. **Primary Perturbation**: A primary in-silico perturbation is applied to the selected TF. The 
#' observed expression levels of the TF are adjusted by simulating changes using a perturbation 
#' direction (`perturb_dir`) and generating a new expression matrix.
#' 
#' 2. **Signal Propagation**: The perturbation signal is propagated throughout the TF regulatory 
#' network using the adjacency matrix derived from TF-target relationships. The propagation process 
#' considers the network structure and propagates the signal over `n_iters` iterations.
#' 
#' 3. **Transition Probability Computation**: Transition probabilities between cells are calculated based on 
#' the perturbation, optionally incorporating velocity or graph topology information.
#' 
#' This method leverages co-expression and regulatory network information to study the potential downstream effects 
#' of transcription factor perturbations on gene expression and cellular states.
#' 
#' @import Seurat
#' @export 
TFPerturbation <- function(
    seurat_obj,
    selected_tf,
    pertub_dir,
    perturbation_name,
    graph='RNA_nn',
    n_iters = 1,
    delta_scale = 1,
    corr_sigma=0.05,
    n_threads=4,
    use_velocyto=TRUE,
    use_graph_tp = FALSE,
    depth = 2,
    target_type = "both",
    use_regulons = TRUE,
    layer = 'counts',
    slot = 'counts',
    assay = 'RNA',
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

    # TODO: add checks

    # get modules 
    modules <- GetModules(seurat_obj, wgcna_name)

    # set up the TF network
    tf_regulons <- GetTFRegulons(seurat_obj)
    tfs <- unique(tf_regulons$tf)
    tf_modules <- subset(modules, gene_name %in% tfs & module != 'grey')

    # get the target genes of the selected tf
    cur_regulon <- subset(tf_regulons, tf == selected_tf)
    cur_regulon_genes <- c( unique(cur_regulon$gene))
    
    # get the TF network for the selected TF
    cur_network <- GetTFTargetGenes(
        seurat_obj,
        selected_tfs=selected_tf, 
        depth=depth, 
        target_type=target_type,
        use_regulons=use_regulons,
        wgcna_name=wgcna_name
    )
    cur_regulon <- subset(cur_network, tf == selected_tf & depth == 1)
    cur_regulon_genes <- c( unique(cur_regulon$gene))

    # convert to igraph
    g <- cur_network %>% 
        dplyr::mutate(score = Gain * sign(Cor)) %>%  
        dplyr::select(c(tf, gene, score)) %>% 
        dplyr::rename(source=tf, target=gene, value=score )
    names(g) <- c('source', 'target', 'value')
    g <- igraph::graph_from_data_frame(
        g, 
        directed=TRUE
    )

    # get adjacency matrix from igraph:
    adj <- igraph::as_adjacency_matrix(g, attr='value')
    tfnet_genes <- rownames(adj)

    # order the genes
    adj <- adj[tfnet_genes, tfnet_genes]
    print(dim(adj))

    # transpose the matrix 
    adj <- t(adj)

    # get cell barcodes:
    cells_use <- colnames(seurat_obj)

    ###########################################################################
    # Set up the observed expression matrix
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

    print(paste0('Applying primary in-silico perturbation on  ', selected_tf))
    exp_per <- ApplyPerturbation(
        seurat_obj,
        exp,
        features = selected_tf,
        perturb_dir = perturb_dir,
        cells_use = cells_use,
        slot = slot,
        assay = assay
    )

    ###########################################################################
    # Part 2: apply signal propagation throughout this module 
    ###########################################################################

    print('Applying signal propagation throughout TF Regulatory network...')
    exp_prop <- ApplyPropagation(
        seurat_obj,
        exp[tfnet_genes,cells_use],
        exp_per[tfnet_genes,cells_use],
        network = adj,
        perturb_dir = perturb_dir,
        delta_scale = delta_scale,
        n_iters = n_iters
    )

    # append the expression matrices:
    exp_simulated <- rbind(
        exp_per[!(rownames(exp_per) %in% tfnet_genes),], # genes that aren't in this module
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
        features=cur_regulon_genes,
        graph=graph, 
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





#' ModulePerturbation
#'
#' This function enables in-silico gene expression perturbation analysis using 
#' a co-expression network. It applies primary perturbations to hub genes, propagates the signal 
#' throughout the co-expression network, and computes cell-cell transition probabilities.
#' 
#' @return A Seurat object containing the in-silico perturbation results as a new assay.
#' 
#' @param seurat_obj A Seurat object containing the gene expression and co-expression data.
#' @param mod Name of the co-expression module to perturb.
#' @param perturb_dir A numeric value determining the type of perturbation to apply: 
#' - negative for knock-down, 
#' - positive for knock-in, 
#' - 0 for knock-out.
#' @param perturbation_name A string representing the name of the in-silico perturbation. 
#' This will be stored as a new assay in the Seurat object.
#' @param graph Name of the cell-cell graph in `Graphs(seurat_obj)`, used for transition probability calculations.
#' @param group.by Optional. A string specifying the column in `seurat_obj@meta.data` used for cell grouping.
#' @param group_name Optional. A string or vector specifying the group(s) within `group.by` to use for perturbation.
#' If NULL, perturbation is applied to all cells.
#' @param n_hubs Number of hub genes to perturb in the selected co-expression module. Default is 5.
#' @param n_iters Number of iterations for propagating the perturbation signal through the network. Default is 3.
#' @param delta_scale A numeric scaling factor controlling the influence of the propagated perturbation. Default is 0.2.
#' @param corr_sigma A numeric scaling factor for adjusting the correlation matrix during transition probability calculations. Default is 0.05.
#' @param n_threads Number of threads to use for parallel computation during correlation calculations. Default is 4.
#' @param use_velocyto Logical. If TRUE, leverages velocyto.R functions for transition probabilities. Default is TRUE.
#' @param use_graph_tp Logical. If TRUE, transition probabilities are computed using the cell-cell graph specified in `graph`. Default is FALSE.
#' @param layer Layer in the assay used for the perturbation. Default is 'counts'.
#' @param slot Slot to extract data for aggregation (e.g., "counts", "data", or "scale.data"). Default is 'counts'.
#' @param assay Name of the assay in `seurat_obj` containing the expression data. Default is 'RNA'.
#' @param wgcna_name Optional. Name of the hdWGCNA experiment in `seurat_obj@misc`. If NULL, defaults to the active WGCNA experiment.
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
    expand_module = 0,
    delta_scale = 0.2,
    corr_sigma=0.05,
    n_threads=4,
    use_velocyto=TRUE,
    use_graph_tp = FALSE,
    use_counts_tp = FALSE,
    layer = 'counts',
    slot = 'counts',
    assay = 'RNA',
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

    # define groups based on group.by
    if(is.null(group.by)){
        group.by <- 'fake_group'
        seurat_obj@meta.data[,group.by] <- "all"
        groups <- c("all")
    } else{
        groups <- unique(as.character(seurat_obj@meta.data[,group.by]))
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

    print(head(modules))

    # Add additional genes to the module? This helps for small modules
    if(expand_module > 0){
    
        cur_kme <- paste0('kME_', mod)
        expand_genes <- modules %>% subset(module == 'grey') %>% 
            slice_max(order_by = get(cur_kme), n=expand_module) %>%
            .$gene_name

        module_genes <- c(module_genes, expand_genes)
    }

    # which cells are we selecting to apply the perturbation?
    cells_use <- colnames(seurat_obj)
    # if(is.null(group.by)){
    #   cells_use <- colnames(seurat_obj)
    # } else{
    #     cells_use <- seurat_obj@meta.data %>% 
    #       subset(get(group.by) %in% group_name) %>%
    #       rownames
    # }

    ###########################################################################
    # Set up the observed expression matrix
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

    print('Applying primary in-silico perturbation to hub genes...')
    exp_per <- ApplyPerturbation(
        seurat_obj,
        exp,
        features = hub_genes,
        perturb_dir = perturb_dir,
        cells_use = cells_use,
        group.by = group.by,
       # group_name = group_name,
        layer = layer,
        slot = slot,
        assay = assay
    )

    ###########################################################################
    # Part 2: apply signal propagation throughout this module 
    ###########################################################################

    # get the TOM for the genes in this module
    cur_TOM <- TOM[module_genes, module_genes]

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

    # experimental: apply propagation separately per group
    # 
    # exp_prop <- do.call(cbind, lapply(groups, function(cur_group){
    #     print(cur_group)
    #     cur_cells <- subset(seurat_obj@meta.data, get(group.by) == cur_group) %>% rownames()
    #     ApplyPropagation(
    #         seurat_obj,
    #         exp[module_genes,cur_cells],
    #         exp_per[module_genes,cur_cells],
    #         network = cur_TOM,
    #         perturb_dir = perturb_dir,
    #         delta_scale = delta_scale,
    #         n_iters = n_iters
    #     )
    # }))

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

    # normalize the perturbed data:
    exp_simulated_norm <- log_normalize(
        exp_simulated, colSums(exp)
    )

    seurat_obj <- SetAssayData(
        seurat_obj, 
        assay = perturbation_name, 
        layer = 'data',
        new.data = exp_simulated_norm
    )

    ###########################################################################
    # Part 3: compute transition probabilities
    ###########################################################################

    if(use_counts_tp){
        layer_tp <- 'counts'; slot_tp <- 'counts'
    } else{
        layer_tp <- 'data'; slot_tp <- 'data'
    }

    print('Computing cell-cell transition probabilities based on the perturbation...')
    seurat_obj <- PerturbationTransitions(
        seurat_obj,
        perturbation_name,
        features=module_genes,
        graph=graph, 
        use_velocyto=use_velocyto,
        use_graph_tp = use_graph_tp,
        corr_sigma=corr_sigma,
        n_threads=n_threads,
        layer=layer_tp,
        slot=layer_tp, 
        assay=assay
    )

    # return the Seurat object
    seurat_obj

}

#' CustomPerturbation
#'
#' This function enables in-silico gene expression perturbation analysis for user-selected 
#' genes using a co-expression network. It applies primary perturbations to the selected 
#' genes, propagates the signal throughout their co-expression neighborhood, and computes 
#' cell-cell transition probabilities.
#'
#' @return A Seurat object containing the in-silico perturbation results as a new assay.
#'
#' @param seurat_obj A Seurat object containing the gene expression and co-expression data.
#' @param selected_features A character vector of gene names to perturb (e.g. candidate driver genes).
#' @param perturb_dir A numeric value determining the type of perturbation to apply: 
#' - negative for knock-down, 
#' - positive for knock-in, 
#' - 0 for knock-out.
#' @param perturbation_name A string representing the name of the in-silico perturbation. 
#' This will be stored as a new assay in the Seurat object.
#' @param graph Name of the cell-cell graph in `Graphs(seurat_obj)`, used for transition probability calculations.
#' @param group.by Optional. A string specifying the column in `seurat_obj@meta.data` used for cell grouping.
#' @param group_name Optional. A string or vector specifying the group(s) within `group.by` to use for perturbation.
#' If NULL, perturbation is applied to all cells.
#' @param n_connections Number of co-expressed genes to include alongside the selected features. 
#' If NULL, defaults to the median module size.
#' @param random_connections Logical. If TRUE, selects random non-selected genes instead of the most 
#' strongly co-expressed genes.
#' @param exclude_grey_genes Logical. If TRUE, excludes grey (unassigned) genes from the co-expression network.
#' @param delta_scale A numeric scaling factor controlling the influence of the propagated perturbation. Default is 0.2.
#' @param corr_sigma A numeric scaling factor for adjusting the correlation matrix during transition probability calculations. Default is 0.05.
#' @param n_threads Number of threads to use for parallel computation during correlation calculations. Default is 4.
#' @param use_velocyto Logical. If TRUE, leverages velocyto.R functions for transition probabilities. Default is TRUE.
#' @param use_graph_tp Logical. If TRUE, transition probabilities are computed using the cell-cell graph specified in `graph`. Default is FALSE.
#' @param use_counts_tp Logical. If TRUE, transition probabilities are computed using the raw counts layer. Default is FALSE.
#' @param layer Layer in the assay used for the perturbation. Default is 'counts'.
#' @param slot Slot to extract data for aggregation (e.g., "counts", "data", or "scale.data"). Default is 'counts'.
#' @param assay Name of the assay in `seurat_obj` containing the expression data. Default is 'RNA'.
#' @param wgcna_name Optional. Name of the hdWGCNA experiment in `seurat_obj@misc`. If NULL, defaults to the active WGCNA experiment.
#'
#' @details 
#' Following co-expression network analysis with hdWGCNA, `GenePerturbation` allows us to perform 
#' in-silico perturbation of any set of user-specified genes. This analysis consists of several key steps:
#'
#' 1. **Apply a primary perturbation** to the selected genes (e.g. knock-in, knock-down, or knock-out).
#'
#' 2. **Apply a secondary perturbation** to their co-expressed neighbors using a signal propagation 
#' algorithm. Given the gene-gene co-expression network from hdWGCNA, the perturbation signal 
#' is propagated iteratively using the network adjacency matrix.
#'
#' 3. **Compute cell-cell transition probabilities** to quantify how the perturbation shifts 
#' cellular states in the latent space or graph.
#'
#' The result is stored as a new assay in the Seurat object with normalized simulated expression data 
#' and computed transition probabilities.
#'
#' @import Seurat
#' @export
CustomPerturbation <- function(
    seurat_obj,
    selected_features,
    perturb_dir,
    perturbation_name,
    graph,
    group.by = NULL,
    group_name = NULL,
    n_connections = NULL,
    random_connections = FALSE,
    exclude_grey_genes = FALSE,
    delta_scale = 0.2,
    corr_sigma=0.05,
    n_threads=4,
    use_velocyto=TRUE,
    use_graph_tp = FALSE,
    use_counts_tp = FALSE,
    layer = 'counts',
    slot = 'counts',
    assay = 'RNA',
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

    # define groups based on group.by
    if(is.null(group.by)){
        group.by <- 'fake_group'
        seurat_obj@meta.data[,group.by] <- "all"
        groups <- c("all")
    } else{
        groups <- unique(as.character(seurat_obj@meta.data[,group.by]))
    }

    # checks 
    # TODO: which checks should go here and which should go inside the functions?
    # it would be nice to do them all here so it errors out immediately.

    # get the TOM:
    TOM <- GetTOM(seurat_obj, wgcna_name)

    # get modules 
    modules <- GetModules(seurat_obj, wgcna_name)

    if(exclude_grey_genes){
        modules <- subset(modules, module != 'grey') %>% 
            dplyr::mutate(module = droplevels(module))
        TOM <- TOM[modules$gene_name, modules$gene_name]
    }

    # calculate the size of each module 
    if(is.null(n_connections)){
        mod_sizes <- table(modules$module)
        mod_sizes <- mod_sizes[names(mod_sizes) != 'grey']
        n_connections <- median(mod_sizes)
    }

    # are we selecting specific genes or random genes?
    if(random_connections){
        print(paste0("Randomly selecting", n_connections, " genes ..."))
        hub_genes <- selected_features
        non_hub_genes <- sample(setdiff(rownames(TOM), hub_genes), n_connections)
        module_genes <- c(hub_genes, non_hub_genes)
    } else{
        print("Using specific features")
        # what are the top connected genes to our selected genes?
        if(length(selected_features == 1)){
            cur_order <- rev(order(TOM[selected_features,]))
        } else{

            # TODO: should this be colSums, or should we calculate median / mean???
            # how similar are the gene sets if we select by sum, median, mean?
            cur_order <- rev(order(colSums(TOM[selected_features,])))
        }
        module_genes <- c(selected_features, rownames(TOM)[cur_order][1:n_connections])
        hub_genes <- selected_features 
        non_hub_genes <- setdiff(module_genes, hub_genes)
    }

    # define the groups
    groups <- unique(as.character(seurat_obj@meta.data[,group.by]))

    # which cells are we selecting to apply the perturbation?
    cells_use <- colnames(seurat_obj)
    # if(is.null(group.by)){
    #   cells_use <- colnames(seurat_obj)
    # } else{
    #     cells_use <- seurat_obj@meta.data %>% 
    #       subset(get(group.by) %in% group_name) %>%
    #       rownames
    # }

    ###########################################################################
    # Set up the observed expression matrix
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

    print('Applying primary in-silico perturbation to selected genes...')
    exp_per <- ApplyPerturbation(
        seurat_obj,
        exp,
        features = hub_genes,
        perturb_dir = perturb_dir,
        cells_use = cells_use,
        group.by = group.by,
       # group_name = group_name,
        layer = layer,
        slot = slot,
        assay = assay
    )

    ###########################################################################
    # Part 2: apply signal propagation throughout this module 
    ###########################################################################

    # get the TOM for the genes in this module
    cur_TOM <- TOM[module_genes, module_genes]

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

    # experimental: apply propagation separately per group
    # 
    # exp_prop <- do.call(cbind, lapply(groups, function(cur_group){
    #     print(cur_group)
    #     cur_cells <- subset(seurat_obj@meta.data, get(group.by) == cur_group) %>% rownames()
    #     ApplyPropagation(
    #         seurat_obj,
    #         exp[module_genes,cur_cells],
    #         exp_per[module_genes,cur_cells],
    #         network = cur_TOM,
    #         perturb_dir = perturb_dir,
    #         delta_scale = delta_scale,
    #         n_iters = n_iters
    #     )
    # }))

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

    # normalize the perturbed data:
    exp_simulated_norm <- log_normalize(
        exp_simulated, colSums(exp)
    )

    seurat_obj <- SetAssayData(
        seurat_obj, 
        assay = perturbation_name, 
        layer = 'data',
        new.data = exp_simulated_norm
    )

    ###########################################################################
    # Part 3: compute transition probabilities
    ###########################################################################

    if(use_counts_tp){
        layer_tp <- 'counts'; slot_tp <- 'counts'
    } else{
        layer_tp <- 'data'; slot_tp <- 'data'
    }

    print('Computing cell-cell transition probabilities based on the perturbation...')
    seurat_obj <- PerturbationTransitions(
        seurat_obj,
        perturbation_name,
        features=module_genes,
        graph=graph, 
        use_velocyto=use_velocyto,
        use_graph_tp = use_graph_tp,
        corr_sigma=corr_sigma,
        n_threads=n_threads,
        layer=layer_tp,
        slot=layer_tp, 
        assay=assay
    )

    # return the Seurat object
    seurat_obj

}

#' @keywords internal
#' @noRd
#'
#' @title log_normalize
#'
#' @description 
#' Internal helper function for normalizing and log-transforming gene expression matrices. 
#' Automatically supports both dense and sparse matrices (`dgCMatrix`).
#'
#' @param X A numeric matrix or sparse matrix of gene expression values (genes x cells).
#' @param col_sums A numeric vector of total counts per cell (length equals number of columns in `X`).
#' @param scale.factor A numeric scaling factor applied after normalization (default = 1e4).
#'
#' @return A numeric or sparse matrix of normalized and log-transformed gene expression values.
#'
#' @details 
#' - For dense matrices, normalization is done via `sweep()`.  
#' - For sparse matrices, normalization is done via efficient column scaling using `Matrix` operations.  
#' - The function applies `log1p()` (log(x + 1)) transformation to all values.  
#'
#' 
log_normalize <- function(X, col_sums, scale.factor = 1e4) {
  if (inherits(X, "dgCMatrix")) {
    # Efficiently normalize sparse matrix columns
    inv_col_sums <- scale.factor / col_sums
    X_norm <- X
    X_norm@x <- X_norm@x * rep(inv_col_sums, diff(X_norm@p))
    # Apply log1p to non-zero entries only
    X_norm@x <- log1p(X_norm@x)
    return(X_norm)
  } else {
    # Dense matrix normalization
    X_norm <- sweep(X, 2, col_sums, "/") * scale.factor
    X_log <- log1p(X_norm)
    return(X_log)
  }
}







# log_normalize <- function(X, col_sums, scale.factor = 1e4) {

#   # Normalize each gene expression value by the total UMI count for its cell
#   X_norm <- sweep(X, 2, col_sums, "/") * scale.factor

#   # Log-transform the normalized values
#   X_log <- log1p(X_norm)
  
#   return(X_log)
# }
