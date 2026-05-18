# FindConnectedComponents

Inspects an SNN (or any named) graph from a Seurat object and labels
each cell by its connected component. The result is written to a new
column in `seurat_obj@meta.data`. This is a useful diagnostic before
running `PerturbationTransitions`, because cells that belong to
different connected components cannot exchange transition probability
mass through the graph.

## Usage

``` r
FindConnectedComponents(
  seurat_obj,
  graph = NULL,
  meta_data_name = "connected_component",
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing at least one graph in `seurat_obj@graphs`.

- graph:

  Character. Name of the graph in `Graphs(seurat_obj)` to analyse. If
  `NULL` (default), the function auto-detects a graph whose name ends in
  `"_snn"`. When multiple SNN graphs exist the first one is used and a
  message is emitted. If no SNN graph is present the first available
  graph is used instead (with a warning).

- meta_data_name:

  Character. Name of the new column written to `seurat_obj@meta.data`.
  Default `"connected_component"`.

- verbose:

  Logical. Print a summary of the component structure. Default `TRUE`.

## Value

The Seurat object with a new integer-factor column in
`seurat_obj@meta.data`. Levels are ordered by component size (largest
component = level 1) so the dominant component always has the lowest
label.

## Details

The graph is treated as undirected for component finding — any non-zero
edge weight is interpreted as a connection. The underlying computation
uses
[`igraph::components()`](https://r.igraph.org/reference/components.html),
which implements a fast depth-first-search.
