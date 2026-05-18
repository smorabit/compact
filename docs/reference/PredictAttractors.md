# Predict Perturbation Attractor States

Performs unsupervised discovery of terminal attractor states (sinks)
within a simulated perturbation vector field. Mathematically, this
function calculates the stationary distribution of the directed Markov
transition matrix by finding the dominant left eigenvector. Cells with
the highest stationary probabilities represent deep phenotypic basins
where the perturbation trajectory naturally pools .

## Usage

``` r
PredictAttractors(
  seurat_obj,
  perturbation_name,
  graph,
  output_name = "attractor_score",
  quantile_threshold = 0.98,
  return_seurat = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing the required graphs.

- perturbation_name:

  Character. The base name of the simulated perturbation. The function
  expects to find a transition probability graph named
  `paste0(perturbation_name, '_tp')` in `seurat_obj@graphs`.

- graph:

  Character. The name of the K-Nearest Neighbors (KNN) graph stored in
  `seurat_obj@graphs` (e.g., "RNA_nn") used to mask the transition
  probabilities.

- output_name:

  Character. The column name used to store the resulting attractor
  scores in `seurat_obj@meta.data`. Default is `"attractor_score"`.

- quantile_threshold:

  Numeric. The percentile threshold used to explicitly define which
  cells belong to the terminal sink state. Default is `0.98` (top 2%).

- return_seurat:

  Logical. If `TRUE`, returns the updated Seurat object. If `FALSE`,
  returns a list containing the numeric attractor scores and a character
  vector of the identified sink cells. Default is `TRUE`.

## Value

If `return_seurat = TRUE`, a Seurat object with the calculated attractor
scores appended to the metadata. If `return_seurat = FALSE`, a list
containing `attractor_score` and `sink_cells`.
