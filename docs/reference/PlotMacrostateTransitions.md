# PlotMacrostateTransitions

Visualizes the coarse-grained group-level transition probability matrix
computed by
[`MacrostateTransitions`](https://smorabit.github.io/compact/reference/MacrostateTransitions.md)
as a heatmap, analogous to the macrostate transition probability heatmap
in CellRank (Lange et al. 2022, Fig. 2c).

Each tile `[source, dest]` is colored by the average transition
probability flowing from the source group to the destination group under
the specified perturbation. The diagonal tiles represent
self-transitions (stability index) and are optionally annotated with the
numeric stability value and a heavier border for visual emphasis.

## Usage

``` r
PlotMacrostateTransitions(
  seurat_obj,
  perturbation_name,
  result_name = NULL,
  group_order = NULL,
  color_scale = c("white", "#2166AC"),
  show_stability = TRUE,
  stability_size = 3,
  diagonal_border = TRUE,
  text_size = 10,
  legend_title = "Transition\nProbability",
  title = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object with results stored in
  `seurat_obj@misc$MacrostateTransitions` (produced by
  [`MacrostateTransitions`](https://smorabit.github.io/compact/reference/MacrostateTransitions.md)).

- perturbation_name:

  Character. Name of the perturbation to visualize. Used to look up
  `result_name` in `@misc$MacrostateTransitions`.

- result_name:

  Character. Key to retrieve from
  `seurat_obj@misc$MacrostateTransitions`. Defaults to
  `perturbation_name`.

- group_order:

  Character vector. Order of groups on both axes. Defaults to
  alphabetical order (the order stored in the result). Both the x-axis
  (destination) and the y-axis (source, reversed) follow this order.

- color_scale:

  Character vector. Two or more colors defining the low-to-high color
  gradient, passed to
  [`ggplot2::scale_fill_gradientn`](https://ggplot2.tidyverse.org/reference/scale_gradient.html).
  Default: `c("white", "#2166AC")` (white to dark blue). For a diverging
  scale (e.g., when plotting ΔQ between perturbations), supply three
  colors with the midpoint as the second element.

- show_stability:

  Logical. If `TRUE` (default), overlay the numeric stability index on
  each diagonal tile.

- stability_size:

  Numeric. Font size for the stability text annotation. Default: `3`.

- diagonal_border:

  Logical. If `TRUE` (default), draw a heavier black border around the
  diagonal tiles to visually highlight self-transitions.

- text_size:

  Numeric. Base text size for axis labels and theme. Default: `10`.

- legend_title:

  Character. Title for the fill legend. Default:
  `"Transition\nProbability"`.

- title:

  Character. Plot title. Defaults to `perturbation_name`.

## Value

A `ggplot2` object.

## References

Lange, M. et al. CellRank for directed single-cell fate mapping. *Nature
Methods* **19**, 159–170 (2022).
[doi:10.1038/s41592-021-01346-6](https://doi.org/10.1038/s41592-021-01346-6)

## See also

[`MacrostateTransitions`](https://smorabit.github.io/compact/reference/MacrostateTransitions.md),
[`PredictAttractors`](https://smorabit.github.io/compact/reference/PredictAttractors.md),
[`VectorFieldCoherence`](https://smorabit.github.io/compact/reference/VectorFieldCoherence.md)
