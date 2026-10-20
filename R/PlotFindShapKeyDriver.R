
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
#' @param out_dir Output directory for the saved plot. Created if it doesn't exist (default: NULL).
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
#' \dontrun{
#' # Requires a Seurat object processed with FindShapKeyDriver().
#' BarplotShap(seurat_obj, top_n = 30)
#'
#' # Save as a PDF with a custom color
#' BarplotShap(
#'   seurat_obj,
#'   save_plot = TRUE,
#'   filename = "barplot_top_shap_DriverGenes.pdf",
#'   bar_color = "#F96815"
#' )
#' }

BarplotShap <- function(seurat_obj,
                        top_n = 20,
                        bar_color = "steelblue",
                        save_plot = FALSE,
                        out_dir = NULL,
                        filename = "barplot_top_shap_DriverGenes.pdf",
                        width = 8,
                        height = 6,
                        dpi = 300) {
#
  # -----------------------------
  # 1) Validate inputs / prerequisites
  # -----------------------------
  # We expect FindShapKeyDriver() to have written shap_summary to seurat_obj@misc$shap.
  # If it isn't present, fail early with a clear error (call.=FALSE keeps errors clean in packages).

  if (is.null(seurat_obj@misc$shap$shap_summary)) {
    stop("shap_summary not found in seurat_obj@misc$shap. Did you run FindShapKeyDriver()?",
         call. = FALSE)
  }

  # Convert to a plain data.frame to avoid S3 class quirks (data.table/tibble/etc.)
  # This makes downstream ggplot evaluation more consistent in package code.
  plot_data <- as.data.frame(seurat_obj@misc$shap$shap_summary)

  # Sanity check: must have at least one row (one gene).
  n_genes <- nrow(plot_data)
  if (n_genes == 0) stop("shap_summary is empty. No genes to plot.", call. = FALSE)

  # -----------------------------
  # 2) Normalize top_n
  # -----------------------------
  if (top_n > n_genes) top_n <- n_genes
  if (top_n < 1) top_n <- min(20, n_genes)

  # Take the first top_n rows. (Assumes shap_summary is already sorted decreasing by importance.)
  plot_data <- utils::head(plot_data, top_n)

  # -----------------------------
  # 3) Choose which SHAP column to plot
  # -----------------------------
  shap_col <- if ("mean_abs_shap" %in% names(plot_data)) "mean_abs_shap" else "mean_shap"

  # -----------------------------
  # 4) Build ggplot (package-safe NSE)
  # -----------------------------
  # Using rlang::.data avoids non-standard evaluation problems in packages.
  # .data[[shap_col]] lets us map a column by string name without get() and without aes_string().
  # stats::reorder orders genes by their SHAP summary so the barplot is ranked.
  # robust ggplot mapping (no get())
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = stats::reorder(.data$variable, .data[[shap_col]]),
      y = .data[[shap_col]]
    )
  ) +
    ggplot2::geom_bar(stat = "identity", fill = bar_color) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Top", top_n, "Genes by Mean SHAP Value"),
      x = "Gene", y = "Mean |SHAP|"
    ) +
    ggplot2::theme_minimal(base_size = 14)


  # -----------------------------
  # 5) Optional saving to disk
  # -----------------------------
  # -----------------------------
  # 5) Optional saving to disk
  # -----------------------------
  if (isTRUE(save_plot)) {

    # Validate required parameters when saving is requested.
    if (is.null(out_dir) || is.null(filename)) {
      stop("When save_plot = TRUE, provide out_dir and filename.", call. = FALSE)
    }

    # Create output directory if missing.
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    file_path <- file.path(out_dir, filename)
    ext <- tools::file_ext(file_path)

    # For PDF, prefer cairo device for better text embedding and vector output.
    # Use grDevices::cairo_pdf (namespaced) to avoid relying on objects in the search path.
    if (tolower(ext) == "pdf") {
      ggplot2::ggsave(
        file_path,
        plot = p,
        width = width,
        height = height,
        device = grDevices::cairo_pdf
      )
    } else {
      # For raster formats (png, jpg), use dpi.
      ggplot2::ggsave(
        file_path,
        plot = p,
        width = width,
        height = height,
        dpi = dpi
      )
    }
  }

  # Return the ggplot object for interactive use.
  return(p)

}

# BarplotShap <- function(seurat_obj,
#                         top_n = 20,
#                         bar_color = "steelblue",
#                         save_plot = FALSE,
#                         out_dir = NULL,
#                         filename = "barplot_top_shap_DriverGenes.pdf",
#                         width = 8,
#                         height = 6,
#                         dpi = 300) {
#   library(ggplot2)
#   library(data.table)
#
#   # Check shap_summary presence
#   if (is.null(seurat_obj@misc$shap$shap_summary)) {
#     stop("shap_summary not found in seurat_obj@misc$shap. Did you run FindShapKeyDriver()?")
#   }
#   shap_summary <- seurat_obj@misc$shap$shap_summary
#
#   # Validate top_n
#   # Sanity check
#   n_genes <- nrow(shap_summary)
#   if (n_genes == 0) stop("shap_summary is empty. No genes to plot.")
#
#   # Adjust top_n
#   if (top_n > n_genes) {
#     message(sprintf("Only %d genes available. Plotting all of them.", n_genes))
#     top_n <- n_genes
#   } else if (top_n < 1) {
#     warning("top_n must be >= 1. Defaulting to 20.")
#     top_n <- min(20, n_genes)
#   }
#
#   # # Subset for plotting
#   # plot_data <- shap_summary[1:top_n]
#   # NOTE shap_summary is a data.table — [1:top_n] works on data.table, but if it's being treated as a data.frame somewhere, it breaks. The real issue is likely that shap_summary came back as a data.frame rather than a data.table (e.g., after being stored in @misc$shap and retrieved).
#   plot_data <- head(as.data.frame(shap_summary), top_n)
#
#   shap_col <- if ("mean_abs_shap" %in% names(plot_data)) "mean_abs_shap" else "mean_shap"
#
#   # Build ggplot
#   p <- ggplot(plot_data, aes(x = reorder(variable, get(shap_col)), y = get(shap_col))) +
#     geom_bar(stat = "identity", fill = bar_color) +
#     coord_flip() +
#     labs(title = paste("Top", top_n, "Genes by Mean SHAP Value"),
#          x = "Gene", y = "Mean |SHAP|") +
#     theme_minimal(base_size = 14)
#
#   # Save plot if requested
#   if (save_plot) {
#     # Check required parameters
#     if (is.null(out_dir) || is.null(filename) || is.null(width) || is.null(height)) {
#       stop("When save_plot = TRUE, you must provide `out_dir`, `filename`, `width`, and `height`.")
#     }
#
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
#     file_path <- file.path(out_dir, filename)
#     ext <- tools::file_ext(file_path)
#
#     if (tolower(ext) == "pdf") {
#       ggsave(file_path, plot = p, width = width, height = height, device = cairo_pdf)
#     } else {
#       ggsave(file_path, plot = p, width = width, height = height, dpi = dpi)
#     }
#   }
#
#   return(p)
# }







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
  # library(ggplot2)
  # library(data.table)
  # library(SHAPforxgboost)
  # library(ggrastr)

  # Extract shap_long from Seurat object
  if (is.null(seurat_obj@misc$shap$shap_long)) {
    stop("shap_long not found in seurat_obj@misc$shap. Did you run FindShapKeyDriver()?", call. = FALSE)
  }
  shap_long <- seurat_obj@misc$shap$shap_long

  # # Ensure output directory exists
  # if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # # Compute top_n variables by mean absolute SHAP
  # top_genes <- shap_long[, list(mean_shap = mean(abs(value))), by = variable][
  #   order(-mean_shap)
  # ]$variable

  # # Compute top_n variables by mean absolute SHAP
  # dt <- data.table::as.data.table(shap_long)
  # genes <- unique(dt$variable)
  # mean_abs <- vapply(genes, function(g) {
  #   mean(abs(dt$value[dt$variable == g]))
  # }, numeric(1))
  # top_genes <- genes[order(-mean_abs)]

  # # Compute top_n variables by mean absolute SHAP
  # dt <- data.table::as.data.table(shap_long)
  # data.table::setDT(dt)
  # summ <- dt[, list(mean_abs_shap = mean(abs(value))), by = "variable"]
  # summ <- summ[order(-mean_abs_shap)]
  # top_genes <- summ$variable
  dt <- data.table::as.data.table(shap_long)
  data.table::setDT(dt)

  if (!all(c("variable", "value") %in% names(dt))) {
    stop("shap_long must contain columns: variable, value. Found: ",
         paste(names(dt), collapse = ", "), call. = FALSE)
  }

  summ <- dt[, list(mean_abs_shap = mean(abs(get("value")))), by = "variable"]
  data.table::setorder(summ, -mean_abs_shap)
  top_genes <- summ$variable



  n_genes <- length(top_genes)
  if (n_genes == 0) stop("shap_long is empty or malformed.")
  if (top_n > n_genes) {
    message(sprintf("Only %d genes available. Plotting all of them.", n_genes))
    top_n <- n_genes
  } else if (top_n < 1) {
    warning("top_n must be >= 1. Defaulting to 20.")
    top_n <- min(20, n_genes)
  }

  # top_genes <- top_genes[1:top_n]
  top_genes <- top_genes[seq_len(top_n)]

  # # Subset and set factor order
  # shap_top <- shap_long[variable %in% top_genes]
  shap_top <- dt[get("variable") %in% top_genes] # NOTE Changed shap_long[variable %in% top_genes] to dt[get("variable") %in% top_genes] — get() forces standard evaluation of the column name string, avoiding NSE

  # shap_top[, variable := factor(variable, levels = rev(top_genes))]
  shap_top[, variable := factor(variable, levels = top_genes)]

  # # Generate beeswarm plot using SHAPforxgboost
  # p <- shap.plot.summary(shap_top, min_color_bound = min_color, max_color_bound = max_color)
  p <- SHAPforxgboost::shap.plot.summary(shap_top, min_color_bound = min_color, max_color_bound = max_color)

  # # If saving the plot
  # if (save_plot) {
  #   # Validate required parameters
  #   if (is.null(out_dir) || is.null(filename) || is.null(width) || is.null(height)) {
  #     stop("When save_plot = TRUE, please provide `out_dir`, `filename`, `width`, and `height`.")
  #   }
  #
  #   if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  #   file_path <- file.path(out_dir, filename)
  #   ext <- tools::file_ext(file_path)
  #
  #   # Rasterize point layer for PDF output
  #   if (tolower(ext) == "pdf") {
  #     # Replace the point layer with a rasterized version
  #     # The 2nd layer is the point layer
  #     # p$layers[[2]] <- ggrastr::geom_point_rast(
  #     #   mapping = p$layers[[2]]$mapping,
  #     #   # # data = p$layers[[2]]$data,
  #     #   # # stat = "sina",
  #     #   # position = position_dodge(width = 0.7),
  #     #   # size = 0.6,
  #     #   # alpha = 0.6
  #     # )
  #     p$layers[[2]] <- ggrastr::rasterise(p$layers[[2]], dpi = dpi)
  #     ggsave(file_path, plot = p, width = width, height = height, device = cairo_pdf)
  #   } else {
  #     ggsave(file_path, plot = p, width = width, height = height, dpi = dpi)
  #   }
  # }
  #
  # # Also rasterize returned plot (useful even if not saving)
  # p$layers[[2]] <- ggrastr::rasterise(p$layers[[2]], dpi = dpi)

  # # Rasterize returned plot (layer index might change; guard it)
  # if (length(p$layers) >= 2) {
  #   p$layers[[2]] <- ggrastr::rasterise(p$layers[[2]], dpi = dpi)
  # }
  if (length(p$layers) > 0) {
    p$layers <- lapply(p$layers, function(lyr) ggrastr::rasterise(lyr, dpi = dpi))
  }

  if (isTRUE(save_plot)) {
    if (is.null(out_dir) || is.null(filename) || is.null(width) || is.null(height)) {
      stop("When save_plot = TRUE, provide out_dir, filename, width, height.", call. = FALSE)
    }
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    file_path <- file.path(out_dir, filename)
    ext <- tools::file_ext(file_path)

    if (tolower(ext) == "pdf") {
      ggplot2::ggsave(
        file_path, plot = p,
        width = width, height = height,
        device = grDevices::cairo_pdf
      )
    } else {
      ggplot2::ggsave(
        file_path, plot = p,
        width = width, height = height,
        dpi = dpi
      )
    }
  }

  return(p)

}
