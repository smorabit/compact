# 0.1.0 (05-20-2026)
## Added
- `MacrostateTransitions` and `PlotMacrostateTransitions`: coarse-grain the cell-level transition probability matrix into a K × K group-level summary, with a stability index on the diagonal and a heatmap visualization.
- `CombinePerturbAssays`: combine multiple perturbation assays by additive delta superposition into a single composite assay.
- New vignette: **Basics — Simulated Branching Trajectory** (`simulation_tutorial`) covering the full `ModulePerturbation` workflow, graph connectivity checks, expression inspection, vector fields, macrostate transitions, and Markov chain analyses.
- New vignette: **Advanced Use Cases — NSCLC CD8+ T Cells** (`advanced_NSCLC`) demonstrating step-by-step use of `ApplyPerturbation`, `ApplyPropagation`, and `PerturbationTransitions` with a custom gene set, a pseudobulk Pearson correlation network, `CustomPerturbation`, and `CombinePerturbAssays`.
- pkgdown documentation site updated with both new vignettes.

## Changes
- None

# 0.0.3 (02-20-2025)
## Added
- New functions for cell fate analysis: `PredictPerturbationTime`, `PredictAttractors`, `PredictFates`, and `PredictCommitment`.
- `VectorFieldCoherence` function to calculate cell-level local coherence scores for the vector field.

## Changes
- None

# 0.0.2 (09-16-2025)
## Added
- `CustomPerturbation` function.

## Changes
- None

# 0.0.1 (01-11-2024)
## Added
- Initial alpha version of **`compact`**

## Changes
- None