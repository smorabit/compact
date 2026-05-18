# Add log2 fold-change columns between two sets of UCell scores

This internal helper computes log2 fold-change between a set of "pre"
UCell scores and matching "post" scores stored in the Seurat object's
metadata. It searches for columns ending in `pre_suffix` and expects
matching columns where `pre_suffix` is replaced by `post_suffix`. For
each pathway, it adds a new column with suffix `logfc_suffix`:

## Usage

``` r
.add_logfc_scores(
  seurat_obj,
  pre_suffix,
  post_suffix,
  logfc_suffix,
  eps = 1e-10
)
```

## Arguments

- seurat_obj:

  A Seurat object with UCell scores stored in `seurat_obj@meta.data`.

- pre_suffix:

  Character suffix used to identify baseline columns (e.g. `"_M6_RNA"`).

- post_suffix:

  Character suffix used to identify perturbation columns (e.g.
  `"_M6_M6up"`).

- logfc_suffix:

  Character suffix used for the new logFC columns (e.g. `"_M6_logFC"`).

- eps:

  Small numeric constant added to both pre and post scores to avoid
  division by zero in the log-ratio.

## Value

The input Seurat object with additional logFC columns added to
`seurat_obj@meta.data`.

## Details

\$\$logFC = \log_2 \left( \frac{\mathrm{post} +
\varepsilon}{\mathrm{pre} + \varepsilon} \right)\$\$
