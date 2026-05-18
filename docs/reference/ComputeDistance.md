# ComputeDistance

Computes pairwise distances between groups (cell states) in a Seurat
object using one of several distance metrics.

## Usage

``` r
ComputeDistance(
  seurat_object,
  groupby,
  reduction,
  method = c("edist", "euclidean", "spearman"),
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
  `"pca"`, `"harmony"`) used for distance calculations.

- method:

  A character scalar specifying the distance metric to use. Supported
  options are `"euclidean"`, `"edist"` (Energy distance), and
  `"spearman"`.

- verbose:

  Logical; if `TRUE`, prints progress messages.

## Value

A symmetric `data.frame` (groups \\\times\\ groups) containing pairwise
distances between groups.

## Details

Groups are defined by `groupby`. Distances are computed in a specified
low-dimensional representation (`reduction`) that preserves meaningful
transcriptional geometry (e.g., PCA or Harmony). Depending on `method`,
distances are computed between group centroids (Euclidean, Spearman) or
between full cell distributions (Energy distance).

Visualization-oriented nonlinear embeddings (e.g., UMAP or t-SNE) are
not recommended for distance-based calculations.

## Examples

``` r
# Euclidean distance between group centroids
ComputeDistance(seurat_object, groupby = "Diagnosis",
                reduction = "pca", method = "euclidean")
#> Error: object 'seurat_object' not found

# Energy distance between group distributions
ComputeDistance(seurat_object, groupby = "Diagnosis",
                reduction = "pca", method = "edist")
#> Error: object 'seurat_object' not found

# Spearman distance between group centroids
ComputeDistance(seurat_object, groupby = "Diagnosis",
                reduction = "pca", method = "spearman")
#> Error: object 'seurat_object' not found
```
