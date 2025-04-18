% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ApplyPropagation.R
\name{ApplyPropagation}
\alias{ApplyPropagation}
\title{ApplyPropagation}
\usage{
ApplyPropagation(
  seurat_obj,
  exp,
  exp_per,
  network,
  perturb_dir = perturb_dir,
  n_iters = 3,
  delta_scale = 0.2
)
}
\arguments{
\item{seurat_obj}{A Seurat object}

\item{exp}{A features by cells matrix containing the observed expression matrix.}

\item{exp_per}{A features by cells matrix containing the perturbation results from ApplyPerturbation.}

\item{network}{A gene-gene network to apply the signal propagation.}

\item{n_iters}{The number of times to apply the signal propagation.}
}
\value{
A dgCMatrix object containing the updated expression matrix with the applied perturbations
}
\description{
This function applies an in-silico perturbation to selected features in a Seurat object.
}
