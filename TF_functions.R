#' MotifScan
#'
#' Scan gene promoters for transcription factor motifs
#'
#' @return seurat_obj with information about motifs in each gene's promoter 
#'
#' @param seurat_obj A Seurat object
#' @param pfm motif matrix set from JASPAR2020 for example
#' @param EnsDb Ensembl database such as EnsDb.Hsapiens.v86 or EnsDb.Mmusculus.v79
#' @param species_genome species genome name. hg38, mm10, etc...
#' @param wgcna_name name of the WGCNA experiment
#' 
#' @details 
#' MotifScan searches the promoter sequences of protein coding genes for
#' transcription factor (TF) motifs. This function first extracts promoter 
#' information using the ensembldb R package, and then motifmatchr is used 
#' to query these promoters for given sets of motifs based on their position
#' frequency matrices (PFMs). The Seurat object is updated with a matrix of 
#' motif-gene hits, a list of putative target genes for each TF, the list of PFMs,
#' and a table of relevant motif information. 
#' 
#' @import Seurat, Matrix, motifmatchr, GenomeInfoDb, GRanges
#' @export
MotifScan <- function(
    seurat_obj,
    pfm, 
    EnsDb, 
    species_genome, 
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

    # TODO: add checks?
    # check pfm
    # check ensdb
    
    # get a dataframe of just the motif name and the motif ID:
    motif_df <- data.frame(
        motif_name = purrr::map(1:length(pfm), function(i){pfm[[i]]@name}) %>% unlist,
        motif_ID = purrr::map(1:length(pfm), function(i){pfm[[i]]@ID}) %>% unlist
    )

    # get promoter and gene coords:
    gene.promoters <- ensembldb::promoters(EnsDb)
    gene.coords <- ensembldb::genes(EnsDb)

    # subset by protein coding
    gene.promoters <- gene.promoters[gene.promoters$tx_biotype == 'protein_coding']
    gene.coords <- gene.coords[gene.coords$gene_biotype == 'protein_coding']

    # subset by main chromosomes
    gene.promoters <- gene.promoters[as.character(GenomeInfoDb::seqnames(gene.promoters)) %in% c(1:100, 'X','Y')]
    gene.coords <- gene.coords[as.character(GenomeInfoDb::seqnames(gene.coords)) %in% c(1:100, 'X', 'Y')]

    # add the gene name to the promoter object
    gene.promoters$symbol <- gene.coords$symbol[base::match(gene.promoters$gene_id, names(gene.coords))]

    # drop unnecessary chromosomes
    gene.promoters <- GenomeInfoDb::keepSeqlevels(
        gene.promoters, 
        value=levels(droplevels(GenomeInfoDb::seqnames(gene.promoters)))
    )

    # rename seqlevels to add 'chr', 
    old_levels <- levels(GenomeInfoDb::seqnames(gene.promoters))
    #new_levels <- ifelse(old_levels %in% c('X', 'Y'), old_levels, paste0('chr', old_levels))
    new_levels <- paste0('chr', old_levels)
    gene.promoters <- GenomeInfoDb::renameSeqlevels(gene.promoters, new_levels)

    # set the genome (not sure if we NEED to do this...)
    GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gene.promoters)) <- species_genome

    # set up promoters object that only has the necessary info for motifmatchr
    my_promoters <- GRanges(
        seqnames =  droplevels(GenomeInfoDb::seqnames(gene.promoters)),
        IRanges(
        start = start(gene.promoters),
        end = end(gene.promoters)
        ),
        symbol = gene.promoters$symbol,
        genome=species_genome
    )

    # scan these promoters for motifs:
    print('Matching motifs...')
    motif_ix <- motifmatchr::matchMotifs(pfm, my_promoters, genome=species_genome)

    # get the matches
    tf_match <- motifmatchr::motifMatches(motif_ix)
    rownames(tf_match) <- my_promoters$symbol

    # use motif names as the column names:
    colnames(tf_match) <- motif_df$motif_name

    # only keep genes that are in the Seurat object and in the given EnsDb:
    gene_list <- rownames(seurat_obj)
    gene_list <- gene_list[gene_list %in% rownames(tf_match)]
    tf_match <- tf_match[gene_list,]

    # get list of target genes for each TF:
    print('Getting putative TF target genes...')
    tfs <- motif_df$motif_name
    tf_targets <- list()
    n_targets <- list()
    for(cur_tf in tfs){
        tf_targets[[cur_tf]] <- names(tf_match[,cur_tf][tf_match[,cur_tf]])
        n_targets[[cur_tf]] <- length(tf_targets[[cur_tf]] )
    }
    n_targets <- unlist(n_targets)

    # add number of target genes to motif_df
    motif_df$n_targets <- n_targets

    # add the gene name to the moti_df
    # remove extra characters from the motif names
    motif_names <- motif_df$motif_name
    tmp <- gsub("\\(.*)", "", motif_names)
    tmp <- gsub('::', ',', as.character(tmp))

    motif_df$tmp <- tmp

    # for motifs that correspond to two genes, split them apart
    tmp <- motif_df$tmp; names(tmp) <- motif_df$motif_ID
    motif_df_tmp <- do.call(rbind,lapply(1:length(tmp), function(i){
    x <- tmp[i]
    id <- names(x)
    if(grepl(',', x)){
        x <- as.character(unlist(do.call(rbind, strsplit(x, ','))))
    }
    data.frame(motif_ID = id, gene_name = as.character(x))
    }))


    # merge with the other motif df:
    ix <- match(motif_df_tmp$motif_ID, motif_df$motif_ID)
    motif_df_tmp <- cbind(motif_df_tmp, motif_df[ix,c('motif_name', 'n_targets')])
    rownames(motif_df_tmp) <- 1:nrow(motif_df_tmp)
    motif_df_tmp <- dplyr::select(motif_df_tmp, c(motif_ID, motif_name, n_targets, gene_name))

    motif_df <- motif_df_tmp

    # subset to only contain genes in the seurat obj
    motif_df <- subset(motif_df, gene_name %in% rownames(seurat_obj))

    # add info to seurat object
    seurat_obj <- SetMotifMatrix(seurat_obj, tf_match)
    seurat_obj <- SetMotifs(seurat_obj, motif_df)
    seurat_obj <- SetMotifTargets(seurat_obj, tf_targets)
    seurat_obj <- SetPFMList(seurat_obj, pfm)

    seurat_obj
}

#' ConstructTFNetwork
#'
#' Construct a network of transcription factors and target genes based on gene co-expression
#'
#' @return seurat_obj with the TF network added
#'
#' @param seurat_obj A Seurat object
#' @param model_params a list of model parameters to pass to xgboost
#' @param nfold number of folds for cross validation
#' @param wgcna_name name of the WGCNA experiment
#' 
#' @details 
#' ConstructTFNetwork uses motif-gene information to build a directed network of transcription 
#' factors (TFs) and target genes. XGBoost regression is leveraged to model the expression of 
#' each gene based on its candidate TF regulators. This analysis gives us information about 
#' how important each TF is for predicting each gene, allowing us to prioritize the most likely
#' regulators of each gene. This process is done on the gene expression matrix stored with SetDatExpr,
#' which is typically the hdWGCNA metacell gene expression matrix.  
#' 
#' @import Seurat, Matrix, xgboost
#' @export
ConstructTFNetwork <- function(
    seurat_obj,
    model_params,
    nfold=5,
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # TODO: add checks? 
   
    # get the motif information from the Seurat object:
    motif_matrix <- GetMotifMatrix(seurat_obj)
    motif_df <- GetMotifs(seurat_obj)

    # check that gene names column is in the motif names
    if(! 'gene_name' %in% colnames(motif_df)){
        stop('gene_name column missing in motif table (GetMotifs(seurat_obj)). Please add a column indicating the gene_name in the seurat_obj for each motif.' )
    }

    # subset the motif_df by genes that are in the seurat obj:
    motif_df <- subset(motif_df, gene_name %in% rownames(seurat_obj))

    # get the expression matrix:
    datExpr <- as.matrix(GetDatExpr(seurat_obj, wgcna_name=wgcna_name))
    genes_use <- colnames(datExpr)

    # set up output dataframes
    importance_df <- data.frame()
    eval_df <- data.frame()

    # set up the progress bar
    pb <- utils::txtProgressBar(min = 0, max = length(genes_use), style = 3, width = 50, char = "=")
    counter <- 1

    ts <- tic()
    for(cur_gene in genes_use){

      print(cur_gene)

        setTxtProgressBar(pb, counter)

        # check if this gene is in the motif matrix:
        if(! cur_gene %in% rownames(motif_matrix)){
            print(paste0('Gene not found in the motif_matrix, skipping ', cur_gene))
            next
        }

        # get the list of TFs that regulate this gene:
        cur_tfs <- names(which(motif_matrix[cur_gene,]))
        cur_tfs <- subset(motif_df, motif_name %in% cur_tfs) %>% .$gene_name %>% unique
        cur_tfs <- cur_tfs[cur_tfs %in% genes_use]

        # set up the expression matrices
        if(cur_gene %in% cur_tfs){
            cur_tfs <- cur_tfs[cur_tfs != cur_gene]
            x_vars <- datExpr[,cur_tfs]
        } 
        x_vars <- datExpr[,cur_tfs]
        y_var <- as.numeric(datExpr[,cur_gene])

        if(length(cur_tfs) < 2){
            print(paste0('Not enough putative TFs, skipping ', cur_gene))
            next
        }

        # correlation:
        tf_cor <- as.numeric(cor(
            x=as.matrix(x_vars),
            y=y_var
        ))
        names(tf_cor) <- cur_tfs

        # run xgboost model
        if(all(y_var == 0)){
            print(paste0('skipping ', cur_gene))
            next
        }
        xgb <- xgboost::xgb.cv(
            params = model_params,
            data = x_vars,
            label = y_var,
            nrounds = 100,
            showsd = FALSE,
            nfold = nfold,
            callbacks = list(cb.cv.predict(save_models=TRUE)),
            verbose=FALSE
        )

        # get the CV evaluation info
        xgb_eval <- as.data.frame(xgb$evaluation_log)
        xgb_eval$variable <- cur_gene

        # average the importance score from each fold
        importance <- Reduce('+', lapply(1:nfold, function(i){
            cur_imp <- xgb.importance(feature_names = colnames(x_vars), model = xgb$models[[i]])
            ix <- match(colnames(x_vars),  as.character(cur_imp$Feature))
            cur_imp <- as.matrix(cur_imp[ix,-1])
            cur_imp[is.na(cur_imp)] <- 0
            cur_imp
        })) / nfold
        importance <- as.data.frame(importance)

        # add tf and source info
        importance$tf <- colnames(x_vars)
        importance$gene<- cur_gene

        # add the tf correlation information
        importance$Cor <- as.numeric(tf_cor)

        # re-order columns, and re-order rows by gain:
        importance <- importance %>% dplyr::select(c(tf, gene, Gain, Cover, Frequency, Cor))
        importance <- arrange(importance, -Gain)

        # append
        importance_df <- rbind(importance_df, importance)
        eval_df <- rbind(eval_df, xgb_eval)

        # update progress bar
        counter <- counter+1

    } 

    # close the progress bar
    te <- toc()
    (te$toc - te$tic) / 60 
    close(pb)

    # add the results to the Seurat object
    seurat_obj <- SetTFNetwork(seurat_obj, importance_df, wgcna_name=wgcna_name)
    seurat_obj <- SetTFEval(seurat_obj, eval_df, wgcna_name=wgcna_name)

    #return(list(importance=importance_df, eval=eval_df))
    seurat_obj

}

#' AssignTFRegulons
#'
#' Define the set of likely target genes (Regulons) for each transcrition factor
#'
#' @return seurat_obj with the TF Regulon information added
#'
#' @param seurat_obj A Seurat object
#' @param strategy method for defining regulons, "A", "B", or "C". See Details for more info.
#' @param reg_thresh threshold for regulatory score in strategies A, B, and C
#' @param n_tfs for strategy A, the number of top TFs to keep for each gene
#' @param n_genes for strategy B, the number of top target genes to keep for each TF
#' @param wgcna_name name of the WGCNA experiment
#' @details 
#' AssignTFRegulons uses the TF network information from ConstructTFNetwork to define 
#' sets of confident TF-gene pairs. A "regulon" is the set of target genes for a given TF,
#' and this function provides three different strategies to define TF regulons. Strategy "A"
#' selects the top TFs for each gene, strategy "B" selects the top genes for each TF, and 
#' strategy "C" retains all TF-gene pairs above a certain regulatory score (reg_thresh). 
#' 
#' @import Seurat, Matrix
#' @export
AssignTFRegulons <- function(
    seurat_obj,
    strategy = "A", # A, B, or C
    reg_thresh = 0.01,
    n_tfs = 10,
    n_genes = 50,
    wgcna_name=NULL
){

    # get data from active assay if wgcna_name is not given
    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # TODO: Checks 

    # get the tf_net from the seurat_obj 
    tf_net <- GetTFNetwork(seurat_obj, wgcna_name)

    if(strategy == 'A'){

        # Take the top TFs for each gene
        tf_regulons <- tf_net %>% 
            subset(Gain >= reg_thresh) %>% 
            group_by(gene) %>%
            slice_max(order_by=Gain, n=n_tfs) %>% 
            ungroup()

    } else if(strategy == 'B'){

        # Take the top target genes for each TF
        tf_regulons <- tf_net %>% 
            subset(Gain >= reg_thresh) %>% 
            group_by(tf) %>%
            slice_max(order_by=Gain, n=n_genes) %>% 
            ungroup()

    } else if(strategy == 'C'){

        # Take all interactions above a certain score
        tf_regulons <- tf_net %>% 
            subset(Gain >= reg_thresh) 

    } else {
        stop('Invalid choice for strategy. Valid choices are A, B, or C.')
    }

    # add the regulons to the seurat object
    seurat_obj <- SetTFRegulons(seurat_obj, tf_regulons, wgcna_name)

}


#' RegulonScores
#'
#' Calculate expression scores for TF Regulons
#'
#' @return seurat_obj with the TF Regulon Scores added
#'
#' @param seurat_obj A Seurat object
#' @param target_type Which types of TF target genes to compute scores for? "positive", "negative", or "both"
#' @param cor_thresh threshold for TF-gene correlation for genes to be included in the regulon score
#' @param exclude_grey_genes option to exclude genes that are in the grey module from the regulon scores
#' @param wgcna_name name of the WGCNA experiment
#' @details 
#' RegulonScores calculates expression signatures for each TF regulon using the UCell algorithm
#' This function can calculate separate scores for TF regulons for target genes that are positively
#' or negatively correlated with the TF, representing putative acivated or repressed genes. These scores 
#' conveniently summarize the expression levels of the entire TF regulon, similar to module eigengenes
#' for the co-expression network analyssis. 
#' 
#' @import Seurat, Matrix, UCell
#' @export
RegulonScores <- function(
    seurat_obj,
    target_type = 'positive', # 'negative', 'all'
    cor_thresh = 0.05,
    exclude_grey_genes = TRUE,
    wgcna_name = NULL,
    ... # options to pass to UCELL
){

    # get modules:
    modules <- GetModules(seurat_obj, wgcna_name)

    # get TF target genes:
    tf_regulons <- GetTFRegulons(seurat_obj, wgcna_name)

    # subset by type of target gene?
    if(target_type == 'positive'){
        tf_regulons <- subset(tf_regulons, Cor > cor_thresh)
    } else if(target_type == 'negative'){
        tf_regulons <- subset(tf_regulons, Cor < cor_thresh)
    } else if(target_type == 'both'){
        tf_regulons <- subset(tf_regulons, abs(Cor) > cor_thresh)
    }

    if(exclude_grey_genes){

        # subset modules
        modules <- subset(modules, module != 'grey')

        # subset regulons
        tf_regulons <- subset(tf_regulons, gene %in% modules$gene_name & tf %in% modules$gene_name)

    }

    # set up the lists
    tfs_use <- unique(tf_regulons$tf)
    target_genes <- lapply(tfs_use, function(cur_tf){
        subset(tf_regulons, tf == cur_tf) %>% .$gene
    })
    names(target_genes) <- tfs_use

    # use UCell to comptue the TF regulons cores
    regulon_scores <- UCell::AddModuleScore_UCell(
        seurat_obj, features=target_genes,
        ...
    )@meta.data
    regulon_scores <- regulon_scores[,paste0(tfs_use, '_UCell')]

    # rename the columns to remove "_UCell"
    colnames(regulon_scores) <- gsub("_UCell", "", colnames(regulon_scores))

    # add the regulon scores to the seurat object
    seurat_obj <- SetRegulonScores(
        seurat_obj, 
        regulon_scores,
        target_type,
        wgcna_name
    )

    seurat_obj
}


#' SetRegulonScores
#'
#' @param seurat_obj A Seurat object
#' @param regulon_scores dataframe storing the TF regulon scores 
#' @param target_type dataframe storing the TF regulon scores 
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetRegulonScores <- function(seurat_obj, regulon_scores, target_type, wgcna_name=NULL){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)


    # if regulon scores have not been set, make a list to store them
    if(is.null(seurat_obj@misc[[wgcna_name]]$regulon_scores)){
        tmp <- list(regulon_scores); names(tmp) <- target_type
        seurat_obj@misc[[wgcna_name]]$regulon_scores <- tmp
    } else{
        seurat_obj@misc[[wgcna_name]]$regulon_scores[[target_type]] <- regulon_scores
    }
    seurat_obj
}

#' GetRegulonScores
#'
#' @param seurat_obj A Seurat object
#' @param target_type dataframe storing the TF regulon scores 
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetRegulonScores <- function(seurat_obj, target_type, wgcna_name=NULL){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)
    
    # get the regulon scores
    seurat_obj@misc[[wgcna_name]]$regulon_scores[[target_type]] 
}


############################
# getters and setters
###########################

#' SetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param tf_net dataframe storing the TF network info in ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFNetwork <- function(seurat_obj, tf_net, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_net <- tf_net
  seurat_obj
}

#' GetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFNetwork <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_net 
}


#' SetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param tf_eval dataframe storing the TF network evaluation info from ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFEval <- function(seurat_obj, tf_eval, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_eval <- tf_eval
  seurat_obj
}

#' GetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFEval <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_eval
}

#' SetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param tf_regulons dataframe storing the TF regulon info from AssignTFRegulons
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFRegulons <- function(seurat_obj, tf_regulons, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_regulons <- tf_regulons
  seurat_obj
}

#' GetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFRegulons <- function(seurat_obj, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_regulons

}

#' GetAvgModuleExpr
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetAvgModuleExpr <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$avg_modules
}


