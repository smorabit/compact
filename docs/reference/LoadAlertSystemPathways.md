# Load AlertSystem pathway gene sets (local or GitHub)

This loads pathway gene sets from a local directory (if `path` is
provided) or automatically downloads them from GitHub if `path = NULL`.

## Usage

``` r
LoadAlertSystemPathways(
  path = NULL,
  genetype = c("mouse", "human"),
  database = "MSigDBhallmark"
)
```

## Arguments

- path:

  Optional local directory OR full file path. If `NULL`, the function
  downloads the gene sets automatically from GitHub.

- genetype:

  Either `"mouse"` or `"human"`.

- database:

  Prefix used in filenames. Default: `"MSigDBhallmark"`.

## Value

A named list of gene sets, where each element is a character vector of
gene symbols corresponding to a pathway.

## Details

Filenames follow the pattern: `<database>_gene_sets_<genetype>.rds`

## Examples

``` r
if (FALSE) { # \dontrun{
# Load pathway gene sets (local directory or GitHub fallback)
alertdatabasePATH <- "/dataPATH/"
genetype <- "mouse"  # or "human"

# Local
gene_sets_mouse <- LoadAlertSystemPathways(
  path = alertdatabasePATH,
  genetype= "mouse",
  database = "MSigDBhallmark"
)

# GitHub
# Load mouse MSigDB Hallmark pathways (GitHub fallback by default)
gene_sets_mouse <- LoadAlertSystemPathways(
  genetype = "mouse",
  database = "MSigDBhallmark"
)

# Load human MSigDB Hallmark pathways
gene_sets_human <- LoadAlertSystemPathways(
  genetype = "human",
  database = "MSigDBhallmark"
)

} # }
```
