# Apply In-Silico Perturbation

This function applies an in-silico perturbation (knock-out, knock-down,
or knock-in) to selected features in a Seurat object.

## Usage

``` r
ApplyPerturbation(
  seurat_obj,
  exp,
  features,
  perturb_dir,
  perturb_mode = "zinb",
  cells_use = NULL,
  group.by = NULL,
  layer = "counts",
  slot = "counts",
  assay = "RNA",
  n_workers = 1
)
```

## Arguments

- seurat_obj:

  A Seurat object containing the dataset.

- exp:

  A features-by-cells matrix (typically a sparse matrix) containing the
  observed expression data.

- features:

  Character vector. The selected features to apply the perturbation on.

- perturb_dir:

  Numeric. Determines the type and magnitude of the perturbation.

  - **ZINB mode** (`perturb_mode = "zinb"`): Negative values for
    knock-down, positive for knock-in, and `0` for knock-out. The
    absolute value scales the ZINB-sampled counts before they are added
    to baseline expression.

  - **Multiplicative mode** (`perturb_mode = "multiplicative"`): A
    positive fold-change applied directly to baseline counts (e.g.,
    `2.0` = 2×, `0.5` = half expression). Must be positive and non-zero;
    use `perturb_dir = 0` for knock-out.

- perturb_mode:

  Character. Perturbation model to use. One of:

  - `"zinb"` (default) — models baseline expression with a Zero-Inflated
    Negative Binomial distribution and adds/subtracts ZINB-sampled
    counts. Group-aware; requires `group.by` and ZINB fitting per
    feature.

  - `"multiplicative"` — scales each cell's baseline counts by
    `perturb_dir`. Cell-specific by construction (delta is proportional
    to observed expression), so cells with zero baseline remain at zero.
    Bypasses ZINB fitting entirely — `group.by` and `n_workers` are
    ignored.

- cells_use:

  Character vector. Specific cells to apply the perturbation to. If
  `NULL`, defaults to handling internally.

- group.by:

  Character. Column in `seurat_obj@meta.data` used to group cells for
  modeling. Only used when `perturb_mode = "zinb"`.

- layer:

  Character. The Seurat v5 layer to extract data from. Default is
  `'counts'`.

- slot:

  Character. The Seurat v4 slot to extract data from. Default is
  `'counts'`.

- assay:

  Character. The assay in `seurat_obj` containing expression
  information. Default is `'RNA'`.

- n_workers:

  Integer. Number of parallel workers for fitting and sampling the ZINB
  model across features. Uses fork-based parallelism
  ([`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html)) so
  memory usage does not scale with worker count. Default is `1`
  (serial). Values \> 1 disable the progress bar. Only used when
  `perturb_mode = "zinb"`.

## Value

A `dgCMatrix` object containing the updated expression matrix with the
applied perturbations.

## Details

Two perturbation models are available via `perturb_mode`:

**ZINB mode** (`"zinb"`): Models the baseline expression of each target
feature using a Zero-Inflated Negative Binomial distribution to capture
dropout and overdispersion characteristics typical of scRNA-seq data.
Samples from this distribution, scales by `perturb_dir`, and adds the
result to baseline expression. Because samples are drawn from a
group-level distribution (IID within each group), the resulting delta is
not cell-specific — cells with zero baseline can receive positive counts
during knock-in.

**Multiplicative mode** (`"multiplicative"`): Scales each cell's
baseline counts directly by `perturb_dir` (a fold change). The
perturbation delta is therefore proportional to each cell's observed
expression, making it cell-specific. This is the recommended mode for
knock-in perturbations when downstream signal propagation is used,
because it preserves the correlation structure between up- and
down-regulation experiments: cells that express a gene highly are
affected most, mirroring the behavior of knock-down.

In both modes, any perturbed counts that fall below zero are strictly
bounded to zero, and the result is rounded to maintain integer count
structure.
