# eucldist

Computes pairwise Euclidean distances between group centroids (mean
embedding vectors) in a specified low-dimensional representation of a
Seurat object.

## Usage

``` r
.eucldist(seurat_object, groupby, reduction, verbose = TRUE)
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

A symmetric `data.frame` (groups \\\times\\ groups) of Euclidean
distances between group centroids.

## Details

For each group defined by `groupby`, the centroid is computed as the
column-wise mean of cell embeddings in `reduction`. Euclidean distances
are then computed between all pairs of group centroids.
