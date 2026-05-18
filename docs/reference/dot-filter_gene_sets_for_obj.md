# Filter gene sets to genes present in a Seurat object

This internal helper takes a list of gene sets and restricts each set to
features present in the provided Seurat object. Gene sets with fewer
than `min_genes` genes after filtering are dropped.

## Usage

``` r
.filter_gene_sets_for_obj(gs_list, seu_obj, min_genes = 5)
```

## Arguments

- gs_list:

  Named list of character vectors, each one a gene set.

- seu_obj:

  A Seurat object.

- min_genes:

  Minimum number of genes required to retain a gene set after
  intersecting with `rownames(seu_obj)`.

## Value

A filtered named list of gene sets, containing only genes present in
`seu_obj` and with at least `min_genes` genes.
