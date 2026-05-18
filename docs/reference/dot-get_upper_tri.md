# This function is used internally to extract the upper triangle of a matrix. It replaces the lower triangle of the input matrix with `NA`, which is useful when working with symmetrical matrices such as distance or correlation matrices where only the upper triangle is needed.

This function is used internally to extract the upper triangle of a
matrix. It replaces the lower triangle of the input matrix with `NA`,
which is useful when working with symmetrical matrices such as distance
or correlation matrices where only the upper triangle is needed.

## Usage

``` r
.get_upper_tri(cormat)
```

## Arguments

- cormat:

  A numeric matrix (e.g., a correlation or distance matrix).

## Value

A matrix with `NA` values in the lower triangle.

## Note

This function is for internal use and is not exported.
