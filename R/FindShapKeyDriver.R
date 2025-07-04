# conda activate MPA
required_pkgs <- c("xgboost", "SHAPforxgboost", "ggplot2", "data.table", "ggrastr")
invisible(lapply(required_pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

#' Identify Key Driver Genes Using SHAP and XGBoost
#'
#' This function trains an XGBoost classifier to distinguish a specified case group from others
#' using gene expression data from a Seurat object. SHAP (SHapley Additive exPlanations) values
#' are computed to quantify the contribution of each gene to the model predictions. The top
#' SHAP-ranked genes are returned as key driver candidates.
#'
#' Three model training modes are supported via the `mode` parameter:
#' \itemize{
#'   \item \code{"full"}: Train on the entire dataset (no test set, no cross-validation).
#'   \item \code{"split"}: Randomly split data into training and testing subsets.
#'   \item \code{"cv"}: Perform k-fold cross-validation.
#' }
#' In \code{"split"} mode, test set AUC is computed and stored. In \code{"cv"} mode,
#' per-fold AUCs and mean AUC are computed. In \code{"full"} mode, AUC is not computed to avoid overestimation.
#'
#' SHAP results and metadata are stored in \code{seurat_obj@misc$shap}. Optionally, results are written to disk
#' under \code{out_dir}, with all files documented in an \code{info.txt}.
#'
#' @param seurat_obj A \code{Seurat} object with normalized gene expression and metadata.
#' @param conditions Column name in \code{seurat_obj@meta.data} indicating the condition labels.
#' @param set_case Character. Case label to compare against all other values in \code{conditions}.
#' @param top_n Integer. Number of top SHAP-ranked genes to return (default: 50).
#' @param max_depth Integer. Maximum depth of each XGBoost tree (default: 4).
#' @param eta Numeric. Learning rate (shrinkage) for XGBoost (default: 0.1).
#' @param nrounds Integer. Number of boosting rounds for XGBoost (default: 50).
#' @param nthread Integer. Number of CPU threads to use for XGBoost (default: 2).
#' @param variable_genes Optional character vector of genes to use. If \code{NULL}, uses \code{VariableFeatures(seurat_obj)}.
#' @param out_dir Optional directory to save results. If \code{NULL}, results are not saved to disk.
#' @param mode Character. Training mode: \code{"full"}, \code{"split"}, or \code{"cv"} (default: \code{"full"}).
#' @param train_fraction Numeric. Fraction of cells used for training in \code{"split"} mode (default: 0.8).
#' @param nfold Integer. Number of folds for cross-validation (default: 5).
#' @param seed Integer. Random seed for reproducibility (default: 123).
#'
#' @return A modified \code{Seurat} object with SHAP results stored in \code{@misc$shap}, including:
#' \itemize{
#'   \item \code{model}: Trained XGBoost model (only in \code{"full"} and \code{"split"} modes)
#'   \item \code{shap_result}: Raw SHAP decomposition object (only in \code{"full"} and \code{"split"} modes)
#'   \item \code{shap_long}: Long-format SHAP values (per-cell, per-gene)
#'   \item \code{shap_summary}: Mean absolute SHAP values per gene
#'   \item \code{key_drivers}: Top \code{top_n} SHAP-ranked genes
#'   \item \code{variable_genes}: Genes used in the model
#'   \item \code{test_auc}: AUC on held-out test set (only in \code{"split"} mode)
#'   \item \code{auc_per_fold}: Vector of AUC scores per fold (only in \code{"cv"} mode)
#'   \item \code{mean_auc}: Mean AUC across folds (only in \code{"cv"} mode)
#'   \item \code{mode}: Training mode used
#' }
#'
#' @section Output Files (if \code{out_dir} is provided):
#' \itemize{
#'   \item \code{model.txt}: XGBoost model summary (not in \code{"cv"} mode)
#'   \item \code{shap_result.rds}: SHAP decomposition object (not in \code{"cv"} mode)
#'   \item \code{shap_summary.txt}: Mean absolute SHAP value per gene
#'   \item \code{variable_genes.txt}: List of genes used for model training
#'   \item \code{key_drivers.txt}: Top SHAP-ranked genes
#'   \item \code{shap_long.rds}: Full SHAP values in long format
#'   \item \code{auc.txt}: AUC from held-out test set (only in \code{"split"} mode)
#'   \item \code{auc_per_fold.txt}: AUC per fold (only in \code{"cv"} mode)
#'   \item \code{info.txt}: Summary of saved outputs
#' }
#'
#' @export
#'
#' @examples
#' # Run SHAP analysis in full mode
#' seurat_obj <- FindShapKeyDriver(seurat_obj, conditions = "Diagnosis", set_case = "AD", top_n = 30)
#'
#' # Run in split mode with custom output directory
#' seurat_obj <- FindShapKeyDriver(seurat_obj, conditions = "Diagnosis", set_case = "AD",
#'                                 mode = "split", out_dir = "results/shap/")
#'
#' # Run in 5-fold CV mode
#' seurat_obj <- FindShapKeyDriver(seurat_obj, conditions = "Diagnosis", set_case = "AD",
#'                                 mode = "cv", nfold = 5)

FindShapKeyDriver <- function(seurat_obj, conditions,
                              set_case, top_n = 50,
                              max_depth = 4, eta = 0.1, nrounds = 50,
                              nthread = 2, variable_genes = NULL, out_dir = NULL,
                              mode = "full",
                              train_fraction = 0.8, nfold = 5, seed = 123) {
  require(Seurat)
  require(xgboost)
  require(SHAPforxgboost)
  require(data.table)

  has_pROC <- requireNamespace("pROC", quietly = TRUE)
  if (!has_pROC) {
    warning("pROC not installed: AUC will not be computed.")
  }

  mode <- match.arg(mode, choices = c("full", "split", "cv"))

  if (missing(conditions)) stop("`conditions` is required but was not provided.")
  if (missing(set_case)) stop("`set_case` is required but was not provided.")

  if (!(set_case %in% seurat_obj[[conditions]][, 1])) {
    stop(sprintf("`set_case = '%s'` not found in `%s` column of Seurat object.", set_case, conditions))
  }

  if (is.null(variable_genes)) {
    message("No gene list supplied — using VariableFeatures from Seurat object.")
    variable_genes <- VariableFeatures(seurat_obj)
  } else {
    message("Using user-supplied gene list.")
  }

  all_genes <- rownames(seurat_obj)
  missing_genes <- setdiff(variable_genes, all_genes)
  if (length(missing_genes) > 0) {
    stop(sprintf("%d supplied genes not found in Seurat object: %s",
                 length(missing_genes), paste(head(missing_genes), collapse = ", ")))
  }

  # should add the check right after you define or confirm variable_genes, and
  # before SHAP computation or top_n usage.
  if (length(variable_genes) < top_n) {
    message(sprintf("Only %d variable genes provided. Setting top_n = %d",
                    length(variable_genes), length(variable_genes)))
    top_n <- length(variable_genes)
  }

  message(sprintf("%d genes selected for SHAP analysis.", length(variable_genes)))
  # X <- t(GetAssayData(seurat_obj, slot = "data")[variable_genes, ])
  # NOTE Error Message
  #   Error in t.default(GetAssayData(seurat_obj, slot = "data")[variable_genes,  :
  #   argument is not a matrix
  # This error means that the object you're trying to transpose with t() is not a matrix,
  # likely a sparse matrix or a non-matrix object like a data frame.
  X <- t(as.matrix(GetAssayData(seurat_obj, slot = "data")[variable_genes, ]))

  y <- as.numeric(seurat_obj[[conditions]] == set_case)
  condition_vals <- seurat_obj[[conditions]][, 1]
  rest_cases <- unique(condition_vals[condition_vals != set_case])
  message(sprintf("Creating binary label for condition: '%s' vs rest.", set_case))
  message(sprintf("The 'rest' includes: %s", paste(rest_cases, collapse = ", ")))

  auc <- NA

  if (mode == "cv") {
    set.seed(seed)
    n <- nrow(X)
    fold_ids <- sample(rep(1:nfold, length.out = n))
    all_shap_long <- list()
    all_aucs <- numeric(nfold)

    for (i in 1:nfold) {
      message(sprintf("Running fold %d of %d...", i, nfold))
      test_idx <- which(fold_ids == i)
      train_idx <- setdiff(seq_len(n), test_idx)

      X_train <- X[train_idx, ]
      y_train <- y[train_idx]
      X_test  <- X[test_idx, ]
      y_test  <- y[test_idx]

      dtrain <- xgb.DMatrix(data = X_train, label = y_train)
      dtest  <- xgb.DMatrix(data = X_test)
      params <- list(objective = "binary:logistic", max_depth = max_depth,
                     eta = eta, nthread = nthread)
      model <- xgboost::xgboost(params = params, data = dtrain,
                                nrounds = nrounds, verbose = 0)

      y_pred <- predict(model, dtest)
      if (has_pROC) {
        auc_fold <- pROC::auc(y_test, y_pred)
        all_aucs[i] <- auc_fold
        message(sprintf("Fold %d AUC: %.3f", i, auc_fold))
      }

      shap_result <- SHAPforxgboost::shap.values(model, X_train = X_test)
      shap_long <- SHAPforxgboost::shap.prep(shap_contrib = shap_result$shap_score, X_train = X_test)
      shap_long$fold <- i
      all_shap_long[[i]] <- shap_long
    }

    combined_shap <- data.table::rbindlist(all_shap_long)
    shap_summary <- combined_shap[, .(mean_abs_shap = mean(abs(value))), by = variable]
    shap_summary <- shap_summary[order(-mean_abs_shap)]
    key_drivers <- head(shap_summary$variable, top_n)
    auc <- mean(all_aucs)

    seurat_obj@misc$shap <- list(
      shap_long = combined_shap,
      shap_summary = shap_summary,
      key_drivers = key_drivers,
      variable_genes = variable_genes,
      auc_per_fold = all_aucs,
      mean_auc = auc,
      mode = mode
    )

  } else {
    if (mode == "split") {
      set.seed(seed)
      n <- nrow(X)
      train_idx <- sample(seq_len(n), size = floor(train_fraction * n))
      test_idx <- setdiff(seq_len(n), train_idx)

      X_train <- X[train_idx, ]
      y_train <- y[train_idx]
      X_test  <- X[test_idx, ]
      y_test  <- y[test_idx]

      message(sprintf("Split data into training (%d) and testing (%d) cells", length(train_idx), length(test_idx)))

      dtrain <- xgb.DMatrix(data = X_train, label = y_train)
      dtest <- xgb.DMatrix(data = X_test)
      model <- xgboost(params = list(objective = "binary:logistic",
                                     max_depth = max_depth, eta = eta, nthread = nthread),
                       data = dtrain, nrounds = nrounds, verbose = 0)

      y_pred <- predict(model, dtest)
      if (has_pROC) {
        auc <- pROC::auc(y_test, y_pred)
        message(sprintf("Test AUC: %.3f", auc))
      }

      shap_result <- shap.values(model, X_train = X_test)
      shap_long <- shap.prep(shap_contrib = shap_result$shap_score, X_train = X_test)

    } else {
      dtrain <- xgb.DMatrix(data = X, label = y)
      model <- xgboost(params = list(objective = "binary:logistic",
                                     max_depth = max_depth, eta = eta, nthread = nthread),
                       data = dtrain, nrounds = nrounds, verbose = 0)

      shap_result <- shap.values(model, X_train = X)
      shap_long <- shap.prep(shap_contrib = shap_result$shap_score, X_train = X)
    }

    if (mode != "cv") {
      shap_summary <- data.table(shap_long)[, .(mean_abs_shap = mean(abs(value))), by = variable]
      shap_summary <- shap_summary[order(-mean_abs_shap)]
      key_drivers <- head(shap_summary$variable, top_n)

      seurat_obj@misc$shap <- list(
        model = model,
        shap_result = shap_result,
        shap_long = shap_long,
        shap_summary = shap_summary,
        key_drivers = key_drivers,
        variable_genes = variable_genes,
        test_auc = auc,
        mode = mode
      )
    }
  }

  # Save outputs based on mode
  if (!is.null(out_dir)) {
    shap_dir <- file.path(out_dir, paste0("shap_output_", mode))
    if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)
    message(sprintf("Saving SHAP outputs to: %s", shap_dir))

    info_lines <- c(paste0("SHAP output summary (mode = '", mode, "'):"), "")

    if (mode %in% c("full", "split")) {
      capture.output(str(seurat_obj@misc$shap$model), file = file.path(shap_dir, "model.txt"))

      saveRDS(seurat_obj@misc$shap$model, file = file.path(shap_dir, "model.rds"))

      saveRDS(seurat_obj@misc$shap$shap_result,
              file = file.path(shap_dir, "shap_result.rds"))
      info_lines <- c(info_lines,
                      "- model.txt: XGBoost model structure",
                      "- shap_result.rds: Raw SHAP decomposition object")
    }

    write.table(seurat_obj@misc$shap$shap_summary,
                file = file.path(shap_dir, "shap_summary.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(data.frame(variable_genes = seurat_obj@misc$shap$variable_genes),
                file = file.path(shap_dir, "variable_genes.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(data.frame(key_drivers = seurat_obj@misc$shap$key_drivers),
                file = file.path(shap_dir, "key_drivers.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    saveRDS(seurat_obj@misc$shap$shap_long,
            file = file.path(shap_dir, "shap_long.rds"))

    info_lines <- c(info_lines,
                    "- shap_summary.txt: Mean absolute SHAP score per gene",
                    "- variable_genes.txt: List of genes used as features",
                    "- key_drivers.txt: Top SHAP-ranked genes",
                    "- shap_long.rds: Per-cell SHAP values in long format")

    if (mode == "cv") {
      write.table(data.frame(auc_per_fold = seurat_obj@misc$shap$auc_per_fold),
                  file = file.path(shap_dir, "auc_per_fold.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      info_lines <- c(info_lines,
                      "- auc_per_fold.txt: AUC scores for each fold",
                      "- No model.txt or shap_result.rds saved because multiple models were trained per fold")
    }

    if (mode == "split") {
      writeLines(sprintf("%.5f", seurat_obj@misc$shap$test_auc), file.path(shap_dir, "auc.txt"))
      info_lines <- c(info_lines, "- auc.txt: AUC value from held-out test set")
    }

    writeLines(info_lines, con = file.path(shap_dir, "info.txt"))
    message("SHAP results saved. See info.txt for details.")
  }

  message(sprintf("SHAP driver gene analysis complete (mode = '%s')", mode))
  return(seurat_obj)
}



#
#' Reload SHAP Output from Disk
#'
#' This function reloads previously saved SHAP analysis results from a specified directory
#' (e.g., created by \code{FindShapKeyDriver()}) and attaches them to a Seurat object.
#'
#' It is useful when SHAP results were not stored in \code{@misc} or when working across sessions.
#'
#' @param seurat_obj A \code{Seurat} object to attach SHAP results to.
#' @param shap_dir Path to the directory containing SHAP outputs (e.g., "results/shap_output_split").
#'
#' @return A \code{Seurat} object with SHAP results loaded into \code{seurat_obj@misc$shap}.
#' @export
#'
#' @examples
#' seurat_obj <- ReloadShapOutput(seurat_obj, shap_dir = "results/shap_output_split")

ReloadShapOutput <- function(seurat_obj, shap_dir) {
  if (!dir.exists(shap_dir)) {
    stop(sprintf("Directory not found: %s", shap_dir))
  }

  message(sprintf("Reloading SHAP outputs from: %s", shap_dir))
  shap_list <- list()

  # Load files if present
  if (file.exists(file.path(shap_dir, "shap_summary.txt"))) {
    shap_list$shap_summary <- data.table::fread(file.path(shap_dir, "shap_summary.txt"))
  }

  if (file.exists(file.path(shap_dir, "variable_genes.txt"))) {
    shap_list$variable_genes <- read.table(file.path(shap_dir, "variable_genes.txt"),
                                            header = TRUE, stringsAsFactors = FALSE)$variable_genes
  }

  if (file.exists(file.path(shap_dir, "key_drivers.txt"))) {
    shap_list$key_drivers <- read.table(file.path(shap_dir, "key_drivers.txt"),
                                         header = TRUE, stringsAsFactors = FALSE)$key_drivers
  }

  if (file.exists(file.path(shap_dir, "shap_long.rds"))) {
    shap_list$shap_long <- readRDS(file.path(shap_dir, "shap_long.rds"))
  }

  if (file.exists(file.path(shap_dir, "shap_result.rds"))) {
    shap_list$shap_result <- readRDS(file.path(shap_dir, "shap_result.rds"))
  }

  if (file.exists(file.path(shap_dir, "model.rds"))) {
    shap_list$model <- readRDS(file.path(shap_dir, "model.rds"))
  }

  if (file.exists(file.path(shap_dir, "auc.txt"))) {
    shap_list$test_auc <- as.numeric(readLines(file.path(shap_dir, "auc.txt")))
  }

  if (file.exists(file.path(shap_dir, "auc_per_fold.txt"))) {
    auc_df <- data.table::fread(file.path(shap_dir, "auc_per_fold.txt"))
    shap_list$auc_per_fold <- auc_df$auc_per_fold
    shap_list$mean_auc <- mean(shap_list$auc_per_fold, na.rm = TRUE)
  }

  # Attempt to infer mode from directory name or info.txt
  # Try to extract mode from folder name
  mode_guess <- gsub("^shap_output_?", "", basename(shap_dir))

  # If not a valid mode, fallback to reading from info.txt
  if (!mode_guess %in% c("full", "split", "cv")) {
    info_file <- file.path(shap_dir, "info.txt")
    if (file.exists(info_file)) {
      info_lines <- readLines(info_file)
      mode_line <- grep("mode = '", info_lines, value = TRUE)
      if (length(mode_line) > 0) {
        mode_guess <- sub(".*mode = '([^']+)'.*", "\\1", mode_line[1])
      } else {
        warning("Mode could not be inferred from info.txt.")
      }
    } else {
      warning("info.txt not found. Mode could not be inferred.")
    }
  }

  shap_list$mode <- mode_guess

  seurat_obj@misc$shap <- shap_list
  message("SHAP outputs successfully loaded into seurat_obj@misc$shap.")
  return(seurat_obj)
}
