# edist

Computes pairwise Energy distances between groups (cell states) in a
Seurat object.

## Usage

``` r
.edist(
  seurat_object,
  groupby,
  reduction,
  sample_correct = TRUE,
  squared = FALSE,
  verbose = TRUE
)
```

## Arguments

- seurat_object:

  A `Seurat` object.

- groupby:

  A character scalar giving the name of a metadata column in
  `seurat_object@meta.data` that defines group labels.

- reduction:

  A character scalar giving the name of a dimensional reduction (e.g.,
  `"pca"`, `"harmony"`) used for distance calculations. This should be a
  geometry-preserving representation; visualization-oriented nonlinear
  embeddings (e.g., `"umap"`, `"tsne"`) are not recommended for distance
  computation.

- sample_correct:

  Logical; if `TRUE`, uses the sample-size corrected (U-statistic)
  normalization for within-group terms (dividing by \\N(N-1)\\ rather
  than \\N^2\\). Note that Energy distance can be negative due to
  finite-sample estimation noise even though the population quantity is
  non-negative.

- squared:

  Logical; if `TRUE`, uses squared Euclidean distances in the embedding
  space (matching `metric="sqeuclidean"` in some Python
  implementations). If `FALSE` (default), uses Euclidean distances.

- verbose:

  Logical; if `TRUE`, prints progress messages.

## Value

A symmetric `data.frame` (groups \\\times\\ groups) of pairwise Energy
distances.

## Details

Cells are grouped according to `groupby`. For each pair of groups, the
empirical Energy distance is computed from pairwise distances between
individual cell embeddings in the specified `reduction` space. The
resulting output is a symmetric groups \\\times\\ groups distance
matrix.
