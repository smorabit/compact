## **`compact`**: co-expression module perturbation analysis 

[![R](https://img.shields.io/github/r-package/v/smorabit/COMPACT)](https://github.com/smorabit/COMPACT/tree/main)
[![ISSUES](https://img.shields.io/github/issues/smorabit/COMPACT)](https://github.com/smorabit/COMPACT/issues)

**`compact`** is a framework for performing CO-expression Module Perturbation Analysis in Cellular Transcriptomes. Building off of our previous method [**`hdWGCNA`**](https://smorabit.github.io/hdWGCNA/), **`compact`** applies direct perturbations to co-expression network hub genes, and uses the network structure to propagate the perturbation signal to other linked genes in the network. This framework is highly flexible to perform knock-in, knock-down, or knock-out perturbations on different networks and sets of genes and in different cell lineages, allowing researchers to explore a wide range of strategies mimicking various experimental conditions and interventions.

**`compact`** is under active development, and is currently an **alpha** version. This means that many features of the final package are missing, and current features are subject to change throughout development. The first major release of **`compact`** will coincide with our forthcoming publication in the coming months.

To get started, please install the package and then visit the [co-expression module perturbation analysis tutorial](https://smorabit.github.io/compact/articles/basic_tutorial.html).

## Installation

We recommend creating an R [conda environment](https://docs.conda.io/en/latest/) environment for COMPACT.

```bash
conda create -n COMPACT -c conda-forge r-base r-essentials
conda activate COMPACT
```

Install the required packages in R

```r
# install BiocManager
install.packages("BiocManager")
BiocManager::install()

# install latest version of Seurat from CRAN
install.packages("Seurat")

# install devtools
BiocManager::install("devtools")

# install hdWGCNA, velocyto.R, and compact from devtools
devtools::install_github('smorabit/hdWGCNA', ref='dev')
devtools::install_github("velocyto-team/velocyto.R")
devtools::install_github('smorabit/COMPACT')

```

![](man/figures/COMPACT_Pipeline.png)
