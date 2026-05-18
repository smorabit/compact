# Predict Perturbation Fates (Forward Diffusion)

Simulates the forward diffusion of a specific cell population under a
perturbed vector field. This function models how a starting group of
cells (sources) will transition over time, revealing their terminal fate
bias. It computes the probability distribution of the population across
the manifold after a specified number of Markov steps .

## Usage

``` r
PredictFates(
  seurat_obj,
  perturbation_name,
  graph,
  source_cells = NULL,
  group.by = NULL,
  group_name = NULL,
  output_name = "forward_fate",
  t_steps = 100,
  tolerance = 1e-05,
  rank_transform = TRUE,
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

- source_cells:

  Character vector. An explicit list of cell barcodes to be defined as
  the starting population. If `NULL`, `group.by` and `group_name` must
  be provided.

- group.by:

  Character. The name of a column in `seurat_obj@meta.data` containing
  cluster or group identities. Used to automatically define source
  cells.

- group_name:

  Character vector. The specific identity class(es) within the
  `group.by` column to define as the starting population.

- output_name:

  Character. The column name used to store the resulting fate
  probabilities in `seurat_obj@meta.data`. Default is `"forward_fate"`.

- t_steps:

  Numeric. The maximum number of forward diffusion steps to simulate.
  Default is `100`.

- tolerance:

  Numeric. The convergence tolerance. If the probability distribution
  stops changing between steps (reaches a stationary state early), the
  simulation will stop. Default is `1e-5`.

- rank_transform:

  Logical. If `TRUE`, transforms the raw probabilities into percentile
  ranks (ECDF) for better dynamic range during visualization. Default is
  `TRUE`.

- return_seurat:

  Logical. If `TRUE`, returns the updated Seurat object. If `FALSE`,
  returns a list containing the probabilities and source cells. Default
  is `TRUE`.

## Value

If `return_seurat = TRUE`, a Seurat object with the calculated fate
scores appended to the metadata. If `return_seurat = FALSE`, a list
containing the fate probabilities and the initial source cells.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate where the "Pre-Exhausted" cells go after down-regulating the exhaustion module
seurat_obj <- PredictFates(
  seurat_obj = seurat_obj,
  perturbation_name = "Exhaustion_down",
  graph = "RNA_nn",
  group.by = "seurat_clusters",
  group_name = "7", # Assuming cluster 7 is the pre-exhausted state
  output_name = "fate_from_cluster7",
  t_steps = 100,
  rank_transform = TRUE
)
FeaturePlot(seurat_obj, features = "fate_from_cluster7")
} # }
```
