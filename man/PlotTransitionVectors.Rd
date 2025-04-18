% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotTransitionVectors.R
\name{PlotTransitionVectors}
\alias{PlotTransitionVectors}
\title{Plot Transition Vectors on a Reduced Dimensional Embedding}
\usage{
PlotTransitionVectors(
  seurat_obj,
  perturbation_name,
  color.by,
  reduction = "umap",
  n_threads = 4,
  arrow_scale = 1,
  grid_resolution = 25,
  max_pct = 0.9,
  point_alpha = 0.25,
  point_size = 1
)
}
\arguments{
\item{seurat_obj}{A Seurat object containing single-cell data, including metadata and embeddings.}

\item{perturbation_name}{A character string specifying the name of the perturbation. This is used
to calculate transition vectors based on the perturbed assay.}

\item{color.by}{A character string indicating the metadata column to color the scatter plot points by.}

\item{reduction}{A character string specifying the name of the dimensional reduction to use for the
embedding (default: 'umap').}

\item{n_threads}{An integer specifying the number of threads to use for parallel processing
(default: 4).}

\item{arrow_scale}{A numeric value to scale the transition vectors, adjusting their length on the plot
(default: 1).}

\item{grid_resolution}{An integer specifying the resolution of the grid for overlaying transition
vectors. A higher resolution yields a finer grid, which can provide more detailed vector placement
(default: 25).}

\item{max_pct}{A numeric value between 0 and 1 indicating the maximum percentile of vector lengths
to display. This parameter helps limit the display of outlier vector lengths, making the plot
easier to interpret (default: 0.90).}

\item{point_alpha}{A numeric value controlling the transparency of scatter plot points, with 1 being
fully opaque and 0 being fully transparent (default: 0.25).}

\item{point_size}{A numeric value defining the size of the scatter plot points (default: 1).}
}
\value{
A ggplot object showing the scatter plot of cells colored by the specified metadata and
overlaid with transition vectors indicating cell-state shifts.
}
\description{
This function visualizes cell-state transitions in a Seurat object by plotting transition vectors on a
dimensionality-reduced embedding (e.g., UMAP). It computes transition vectors based on a specified
perturbation and overlays them on a scatter plot of the cells, allowing insights into cell-state shifts
following the perturbation.
}
\details{
This function retrieves the embedding coordinates from the specified dimensional reduction
in the Seurat object. It then calculates the transition vectors using \code{PerturbationVectors}, which
returns vector coordinates (\code{ars}) and distances (\code{arsd}). A grid of arrows is created using
\code{grid_vectors}, allowing for visual simplification of transitions by grouping vectors into a grid
at the specified resolution. The resulting plot includes the scatter plot of cells with color
defined by \code{color.by} and arrows showing the direction and relative magnitude of transitions.
}
