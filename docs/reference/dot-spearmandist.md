# spearmandist

Computes pairwise Spearman distances between group centroids (mean
embedding vectors) in a Seurat object.

## Usage

``` r
.spearmandist(seurat_object, groupby, reduction, verbose = TRUE)
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

- verbose:

  Logical; if `TRUE`, prints progress messages.

## Value

A symmetric `data.frame` (groups \\\times\\ groups) containing pairwise
Spearman distances between group centroids.

## Details

Cells are grouped according to `groupby`. For each group, a centroid is
computed as the column-wise mean of cell embeddings in the specified
`reduction` space. Spearman's rank correlation coefficient is then
computed between all pairs of group centroids and converted to a
distance measure as \\1 - \rho\\.

This metric captures monotonic similarity between group-level expression
profiles and is insensitive to absolute scale, but does not account for
within-group heterogeneity or higher-order distributional differences.
