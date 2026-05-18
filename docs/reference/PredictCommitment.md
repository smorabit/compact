# Predict Fate Commitment (Committor Probabilities)

Calculates the committor probability for every cell under a simulated
perturbation. This function models the vector field as a Dirichlet
boundary value problem to answer a specific question: "What is the exact
probability that a cell will successfully reach the target state (Sink)
before it falls back into the origin state (Source)?" This is highly
effective for identifying epigenetic "points of no return" and
transcriptomic bottlenecks .

## Usage

``` r
PredictCommitment(
  seurat_obj,
  perturbation_name,
  graph,
  source_cells = NULL,
  sink_cells = NULL,
  group.by = NULL,
  source_group = NULL,
  sink_group = NULL,
  output_name = "commitment_score",
  max_iter = 2000,
  tolerance = 1e-05,
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

  Character vector. Cell barcodes to be defined as the starting
  population (Boundary condition = 0).

- sink_cells:

  Character vector. Cell barcodes to be defined as the target population
  (Boundary condition = 1).

- group.by:

  Character. The name of a column in `seurat_obj@meta.data`. Used to
  automatically define source and sink cells if explicit lists are not
  provided.

- source_group:

  Character vector. The identity class(es) in `group.by` representing
  the source state.

- sink_group:

  Character vector. The identity class(es) in `group.by` representing
  the sink state.

- output_name:

  Character. The column name used to store the resulting probabilities
  in `seurat_obj@meta.data`. Default is `"commitment_score"`.

- max_iter:

  Numeric. The maximum number of iterations for the solver. Default is
  `2000`.

- tolerance:

  Numeric. The convergence tolerance. Default is `1e-5`.

- return_seurat:

  Logical. If `TRUE`, returns the updated Seurat object. If `FALSE`,
  returns a list containing the probabilities. Default is `TRUE`.

## Value

If `return_seurat = TRUE`, a Seurat object with the calculated committor
probabilities appended to the metadata. If `return_seurat = FALSE`, a
list containing the probability scores.
