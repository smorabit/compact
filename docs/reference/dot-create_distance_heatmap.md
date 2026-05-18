# This function generates a heatmap from a given matrix. It handles both matrices and data frames (converting them to matrices if needed). It also allows customization of the color palette, axis visibility, and legend display.

This function generates a heatmap from a given matrix. It handles both
matrices and data frames (converting them to matrices if needed). It
also allows customization of the color palette, axis visibility, and
legend display.

## Usage

``` r
.create_distance_heatmap(
  df_matrix,
  title,
  min_val,
  max_val,
  custom_palette = NULL,
  show_x_axis = TRUE,
  show_y_axis = TRUE,
  show_legend = TRUE,
  custom_order = NULL
)
```

## Arguments

- df_matrix:

  A numeric matrix or a data frame to be converted to a matrix.

- title:

  The title for the heatmap.

- min_val:

  The minimum value for the color scale.

- max_val:

  The maximum value for the color scale.

- custom_palette:

  A vector of colors to define the color palette. If `NULL`, a default
  palette is used.

- show_x_axis:

  Logical; whether to display the x-axis labels. Defaults to `TRUE`.

- show_y_axis:

  Logical; whether to display the y-axis labels. Defaults to `TRUE`.

- show_legend:

  Logical; whether to display the legend. Defaults to `TRUE`.

- custom_order:

  A character vector specifying the order of clusters for the x and y
  axes.

## Value

A ggplot2 object representing the heatmap.

## Note

This function is for internal use and is not exported.
