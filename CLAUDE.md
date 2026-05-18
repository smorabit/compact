## Overview
`compact` is an R package to perform in-silico gene expression perturbation analysis in single-cell RNA-seq data, building off our previous method [**`hdWGCNA`**](https://smorabit.github.io/hdWGCNA/). `compact` applies a primary perturbation to a set of selected genes, and leverages the gene co-expression network structure to propagate the signal to other linked genes in the network. `compact` infers the downstream consequences on the cell phenotypic manifold by calculating cell-cell transition probabilities, similar to RNA Velocity, and by performing Markov chain analyses.

## Main funcitonality
The main function of `compact` is `ModulePerturbation`. This function consists of three major steps: `ApplyPerturbation` to simulate the direct perturbation to selected genes, `ApplyPropagation` to propagate the signal to linked genes via the network structure, and `PerturbationTransitions` to calculate the cell-cell transition probabilities. These three steps are highly modular and can be used in different contexts, for instance in the function `TFPertrubation` to perform perturbations of transcription factors and their downstream target genes.  

## Markov chain analyses
`PredictCellFates.R` contains the functions to perform different types of cell-fate analyses.

## Vector field visualizations
Similar to RNA velocity, we can visualize the in-silico perturbation results directly on a two-dimensional view of the scRNA-seq dataset (ie, PCA or UMAP) as a vector field plot (relevant functions in `PlotTransitionVectors.R`)

## Implementation
We implement this as an R package formatted with roxygen2 and pkgdown for the documentation website. `compact` relies on Seurat objects as the underlying data structure.

## Testing
- Tests use testthat via devtools
- Run tests with: devtools::test()
- Check package with: devtools::check()
- Test files live in tests/testthat/
- Prefer testing exported functions over internal helpers
- After writing testing code, ensure to run all tests. If errors are encountered, write solutions to those errors before continuing.
- After testing and debugging, push the changes to Git on the bugfixes branch 

## Code style
- All functions should have detailed docstrings.
- Inline comments should be lowercase (with the exception of acronyms, function names, etc) and without excessive punctuation. 