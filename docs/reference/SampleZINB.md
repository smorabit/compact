# SampleZINB

This function simulated expression values based on a ZINB model that has
been fit to gene expression data.

## Usage

``` r
SampleZINB(model, yobs, ncells = NULL)
```

## Arguments

- model:

  a model object

- ncells:

  number of cells to simulate expression values for. By default will
  simulate data for each cell in the input dataset.

## Value

A numeric vector of simulated expression data based on a ZINB model
