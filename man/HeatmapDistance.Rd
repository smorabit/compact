% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distance_plotting.R
\name{HeatmapDistance}
\alias{HeatmapDistance}
\title{HeatmapDistance
Generate Heatmaps for Original and Perturbed Matrices

This function generates two heatmaps from two matrices (original and perturbed) and displays them side by side on the same color scale.}
\usage{
HeatmapDistance(
  df_original,
  df_perturbed,
  custom_palette = NULL,
  title_original = "Original Assay Cluster Similarity Distance",
  title_perturbed = "Perturbed Assay Cluster Similarity Distance",
  custom_order = NULL
)
}
\arguments{
\item{df_original}{A numeric matrix representing the original (unperturbed) data.}

\item{df_perturbed}{A numeric matrix representing the perturbed data.}

\item{custom_palette}{A vector of colors to define the color palette. Defaults to a red/blue gradient.}

\item{title_original}{The title for the original heatmap. Defaults to "Original Assay Cluster Similarity Distance".}

\item{title_perturbed}{The title for the perturbed heatmap. Defaults to "Perturbed Assay Cluster Similarity Distance".}

\item{custom_order}{A character vector specifying the order of clusters for the x and y axes.}
}
\value{
A patchwork object combining the two heatmaps.
}
\description{
HeatmapDistance
Generate Heatmaps for Original and Perturbed Matrices

This function generates two heatmaps from two matrices (original and perturbed) and displays them side by side on the same color scale.
}
\examples{
p <- HeatmapDistance(df_edist_observed, df_edist_perturbed) # , custom_order = custom_order
}
