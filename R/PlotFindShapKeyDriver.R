# conda activate MPA
required_pkgs <- c("xgboost", "SHAPforxgboost", "ggplot2", "data.table", "ggrastr")
invisible(lapply(required_pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))


#
#' Plot Bar Chart of Top SHAP Driver Genes
#'
#' Generates a horizontal bar plot of the top `n` genes ranked by mean absolute SHAP value,
#' as computed by \code{FindShapKeyDriver()}. The plot provides an interpretable summary of
#' the most important features (genes) driving model predictions.
#'
#' The function pulls SHAP summary statistics from \code{seurat_obj@misc$shap$shap_summary},
#' creates a bar plot of the top `n` genes, and saves it to disk as a PNG or PDF.
#'
#' @param seurat_obj A \code{Seurat} object containing SHAP results from \code{FindShapKeyDriver()},
#'        with a \code{shap_summary} table stored in \code{seurat_obj@misc$shap}.
#' @param top_n Number of top genes to display (default: 20).
#' @param bar_color Fill color for the bars (default: \code{"steelblue"}).
#' @param save_plot Save the Barplot automatically (default: save_plot = FALSE).
#' @param out_dir Output directory for the saved plot. Created if it doesn't exist (default: out_dir = NULL}).
#' @param filename Name of the plot file to save (e.g., \code{"barplot_top_shap_DriverGenes.pdf"}).
#'        File format is determined by the extension (\code{.pdf} or \code{.png}).
#' @param width Plot width in inches (default: 8).
#' @param height Plot height in inches (default: 6).
#' @param dpi Resolution in dots per inch for raster formats like PNG (default: 300).
#'
#' @return A \code{ggplot} object containing the bar plot of SHAP summary values.
#' @export
#'
#' @examples
#' # After running FindShapKeyDriver():
#' BarplotShap(seurat_obj, top_n = 30)
#'
#' # Save as a PDF with a custom color
#' BarplotShap(seurat_obj,save_plot = TRUE,filename = "barplot_top_shap_DriverGenes.pdf", bar_color = "#F96815")

BarplotShap <- function(seurat_obj,
                        top_n = 20,
                        bar_color = "steelblue",
                        save_plot = FALSE,
                        out_dir = NULL,
                        filename = "barplot_top_shap_DriverGenes.pdf",
                        width = 8,
                        height = 6,
                        dpi = 300) {
  library(ggplot2)
  library(data.table)

  # Check shap_summary presence
  if (is.null(seurat_obj@misc$shap$shap_summary)) {
    stop("shap_summary not found in seurat_obj@misc$shap. Did you run FindShapKeyDriver()?")
  }
  shap_summary <- seurat_obj@misc$shap$shap_summary

  # Validate top_n
  # Sanity check
  n_genes <- nrow(shap_summary)
  if (n_genes == 0) stop("shap_summary is empty. No genes to plot.")

  # Adjust top_n
  if (top_n > n_genes) {
    message(sprintf("Only %d genes available. Plotting all of them.", n_genes))
    top_n <- n_genes
  } else if (top_n < 1) {
    warning("top_n must be >= 1. Defaulting to 20.")
    top_n <- min(20, n_genes)
  }

  # Subset for plotting
  plot_data <- shap_summary[1:top_n]
  shap_col <- if ("mean_abs_shap" %in% names(plot_data)) "mean_abs_shap" else "mean_shap"

  # Build ggplot
  p <- ggplot(plot_data, aes(x = reorder(variable, get(shap_col)), y = get(shap_col))) +
    geom_bar(stat = "identity", fill = bar_color) +
    coord_flip() +
    labs(title = paste("Top", top_n, "Genes by Mean SHAP Value"),
         x = "Gene", y = "Mean |SHAP|") +
    theme_minimal(base_size = 14)

  # Save plot if requested
  if (save_plot) {
    # Check required parameters
    if (is.null(out_dir) || is.null(filename) || is.null(width) || is.null(height)) {
      stop("When save_plot = TRUE, you must provide `out_dir`, `filename`, `width`, and `height`.")
    }

    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    file_path <- file.path(out_dir, filename)
    ext <- tools::file_ext(file_path)

    if (tolower(ext) == "pdf") {
      ggsave(file_path, plot = p, width = width, height = height, device = cairo_pdf)
    } else {
      ggsave(file_path, plot = p, width = width, height = height, dpi = dpi)
    }
  }

  return(p)
}







#
#' Plot SHAP Beeswarm Summary for Top Driver Genes
#'
#' Generates a beeswarm-style summary plot of SHAP values for the top `n` driver genes,
#' as computed from `FindShapKeyDriver()`. This visualization helps interpret gene-level
#' contributions to model predictions, showing both SHAP value and feature value per cell.
#'
#' The function extracts SHAP outputs from the `@misc$shap` slot of a Seurat object and
#' plots per-cell SHAP values for the most important genes. The plot is saved as a PDF
#' (with rasterized points) or PNG depending on the filename extension.
#'
#' @param seurat_obj A \code{Seurat} object containing SHAP results from \code{FindShapKeyDriver()},
#'        stored in \code{seurat_obj@misc$shap}.
#' @param top_n Number of top driver genes to display (default: 20).
#' @param save_plot Save the Beeswarmplot automatically (default: save_plot = FALSE).
#' @param out_dir Output directory for the saved plot file. Will be created if it doesn't exist.
#' @param filename Name of the output plot file (e.g., \code{"beeswarm_top_shap_DriverGenes.pdf"}).
#'        File format is inferred from the extension (.pdf or .png).
#' @param min_color Color for low feature values (default: \code{"#FFCC33"}).
#' @param max_color Color for high feature values (default: \code{"#6600CC"}).
#' @param width Width of the plot in inches (default: 8).
#' @param height Height of the plot in inches (default: 6).
#' @param dpi Resolution in dots per inch (only used for raster outputs like PNG; default: 300).
#'
#' @return A \code{ggplot} object representing the beeswarm summary plot.
#' @export
#'
#' @examples
#' # Assuming SHAP analysis was performed using FindShapKeyDriver()
#' BeeswarmplotShap(seurat_obj, out_dir = "figures/", top_n = 30)
#'
#' # To save as a raster-aware PDF
#' BeeswarmplotShap(seurat_obj,save_plot = TRUE,filename = "beeswarm_top_shap_DriverGenes.pdf",out_dir = out_dir)

BeeswarmplotShap <- function(seurat_obj,
                             top_n = 20,
                             save_plot = FALSE,
                             out_dir = NULL,
                             filename = "beeswarm_top_shap_DriverGenes.pdf",
                             min_color = "#FFCC33",
                             max_color = "#6600CC",
                             width = 8,
                             height = 6,
                             dpi = 300) {
  library(ggplot2)
  library(data.table)
  library(SHAPforxgboost)
  library(ggrastr)

  # Extract shap_long from Seurat object
  if (is.null(seurat_obj@misc$shap$shap_long)) {
    stop("shap_long not found in seurat_obj@misc$shap. Did you run FindShapKeyDriver()?")
  }
  shap_long <- seurat_obj@misc$shap$shap_long

  # # Ensure output directory exists
  # if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Compute top_n variables by mean absolute SHAP
  top_genes <- shap_long[, .(mean_shap = mean(abs(value))), by = variable][
    order(-mean_shap)
  ]$variable

  n_genes <- length(top_genes)
  if (n_genes == 0) stop("shap_long is empty or malformed.")
  if (top_n > n_genes) {
    message(sprintf("Only %d genes available. Plotting all of them.", n_genes))
    top_n <- n_genes
  } else if (top_n < 1) {
    warning("top_n must be >= 1. Defaulting to 20.")
    top_n <- min(20, n_genes)
  }

  top_genes <- top_genes[1:top_n]

  # Subset and set factor order
  shap_top <- shap_long[variable %in% top_genes]
  # shap_top[, variable := factor(variable, levels = rev(top_genes))]
  shap_top[, variable := factor(variable, levels = top_genes)]

  # Generate beeswarm plot using SHAPforxgboost
  p <- shap.plot.summary(shap_top, min_color_bound = min_color, max_color_bound = max_color)

  # If saving the plot
  if (save_plot) {
    # Validate required parameters
    if (is.null(out_dir) || is.null(filename) || is.null(width) || is.null(height)) {
      stop("When save_plot = TRUE, please provide `out_dir`, `filename`, `width`, and `height`.")
    }

    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    file_path <- file.path(out_dir, filename)
    ext <- tools::file_ext(file_path)

    # Rasterize point layer for PDF output
    if (tolower(ext) == "pdf") {
      # Replace the point layer with a rasterized version
      # The 2nd layer is the point layer
      # p$layers[[2]] <- ggrastr::geom_point_rast(
      #   mapping = p$layers[[2]]$mapping,
      #   # # data = p$layers[[2]]$data,
      #   # # stat = "sina",
      #   # position = position_dodge(width = 0.7),
      #   # size = 0.6,
      #   # alpha = 0.6
      # )
      p$layers[[2]] <- ggrastr::rasterise(p$layers[[2]], dpi = dpi)
      ggsave(file_path, plot = p, width = width, height = height, device = cairo_pdf)
    } else {
      ggsave(file_path, plot = p, width = width, height = height, dpi = dpi)
    }
  }

  # Also rasterize returned plot (useful even if not saving)
  p$layers[[2]] <- ggrastr::rasterise(p$layers[[2]], dpi = dpi)

  return(p)

}


#
