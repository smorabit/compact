# Predict Perturbation Pseudotime (Absorbing Markov Chain Hitting Time)

Calculates the expected hitting time for cells to reach a defined "sink"
or target state under a simulated perturbation. Mathematically, this
function models the perturbation transition probability matrix as an
Absorbing Markov Chain. It calculates the expected number of random walk
steps required for any transient cell to reach the absorbing sink
states, effectively quantifying the transcriptomic "distance" or
resistance to the perturbation-driven trajectory.

## Usage

``` r
PredictPerturbationTime(
  seurat_obj,
  perturbation_name,
  graph,
  sink_cells = NULL,
  group.by = NULL,
  group_name = NULL,
  output_name = "perturbation_pseudotime",
  max_iter = 10000,
  tolerance = 1e-04,
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
  `seurat_obj@graphs` (e.g., "RNA_nn"). This is used to mask the
  transition probabilities and force the random walk to strictly obey
  the biological manifold.

- sink_cells:

  Character vector. An explicit list of cell barcodes to be defined as
  the target/absorbing states. If `NULL`, `group.by` and `group_name`
  must be provided.

- group.by:

  Character. The name of a column in `seurat_obj@meta.data` containing
  cluster or group identities. Used to automatically define sink cells.

- group_name:

  Character vector. The specific identity class(es) within the
  `group.by` column to define as the sink/target states.

- output_name:

  Character. The column name used to store the resulting pseudotime
  values in `seurat_obj@meta.data`. Default is
  `"perturbation_pseudotime"`.

- max_iter:

  Numeric. The maximum number of iterations for the fixed-point
  iterative sparse solver. Default is `10000`.

- tolerance:

  Numeric. The convergence tolerance for the iterative solver. Default
  is `1e-4`.

- return_seurat:

  Logical. Whether to return a Seurat object (default), or to return the
  calculated perturbation pseudotime as a numeric vector.

## Value

A Seurat object with the calculated perturbation times appended to the
metadata under the column name specified by `output_name`. Sink cells
are strictly assigned a time of 0.

## Details

This function uses a highly efficient, matrix-free iterative solver to
calculate hitting times, ensuring it scales comfortably to hundreds of
thousands of cells without dense matrix inversion. Cells that hit the
`max_iter` ceiling are likely on disconnected manifold "islands" or face
severe topological bottlenecks preventing them from reaching the
specified sink states under the current perturbation.
