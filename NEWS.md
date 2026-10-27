# 0.1.4 (10-27-2026)
## Added
- **Distance-difference heatmaps** for comparing cellular-state-similarity (E-distance) matrices before and after a perturbation. `HeatmapDistanceDiff()` draws a diverging post − pre difference heatmap; `HeatmapDistanceTriptych()` lays out pre, post, and difference panels together (sequential shared scale on pre/post via `HeatmapDistance()`, diverging scale on the difference); and `DiffLimits()` computes one shared colour-limit pair across several difference matrices so they can be read on a common scale.
- New vignette **Perturbation Assessment — Distance-based State Similarity** (`distance_assessment`): quantifying and visualizing how far perturbed cells move from baseline using the distance and distance-difference heatmaps.
- New vignette **SHAP Key Driver Genes** (`shap_key_drivers`): identifying and visualizing the genes that most drive a perturbation's predicted effect with `FindShapKeyDriver()`, `BarplotShap()`, and `BeeswarmplotShap()`.
- A new **Perturbation Assessment** section in the pkgdown navbar grouping the assessment vignettes.

## Changes
- **`CombinePerturbAssays()` reworked for the mixed-layer perturbation workflow.** Direct targets (changed in the `counts` layer) and propagation-driven secondary effects (changed in the log-normalized `data` layer) are now combined differently: count-space deltas are summed onto the baseline counts and re-log-normalized to reproduce the direct component, while each assay's data-space delta is added only for genes it did *not* directly change (component-wise exclusion), so a gene that is a direct target in one assay still keeps its secondary effect from another. Adds provenance checking (`check_provenance`, `require_provenance`, `strict_provenance_check`) and independent count/data flooring (`floor_counts`, `floor_data`).
- **`AlertSystemScore()` substantially reworked.** Two scoring channels are now supported: `expression` (mean post − pre log-normalized expression over retained pathway genes) and `ucell` (post − pre UCell `"difference"`, or the legacy `"log2_ratio"`). Cell-level responses are summarized within `group.by` groups, reporting response prevalence, directionality, and leading genes separately, with optional sample-level sign-flip null calibration (`sample.by`, `n_permutations`). When `out_dir` is supplied, the full state is saved as RDS and the metrics exported as CSV.
- **`AlertHeatmap()` updated** to match the reworked scorer: evidence-aware pathway selection on the positive and negative sides, explicit labeling of UCell difference vs. log2-ratio results (mixed transformations are now rejected rather than silently plotted together), an optional null-standardized (`null_z`) fill, and unit scaling.
- **Alert System save/load refactored** into `R/AlertSystemHelper.R`. `SaveAlertSystem()`, `LoadAlertSystem()`, and `LoadAlertSystemPathways()` now validate the analysis state and gene sets before saving or attaching.
- **`PlotMarkovEmbedding()` / `PlotMarkovScatter()` refined**: cells are layered (unconverged in grey beneath converged), and colour-scale limits are computed from converged cells only, so the scale is not distorted by cells that have not yet arrived.
- `StampPerturbProvenance()` is no longer exported.
- New package dependencies: `dplyr`, `methods`, and `utils` added to `Imports`.
- Expanded installation instructions in the README with separate mamba (recommended) and conda tracks.

# 0.1.3 (08-06-2026)
## Changes
- Fixed `ApplyPropagation` (and therefore `ModulePerturbation`, `TFPerturbation`, `CustomPerturbation`, `TargetModulePerturbation`) diluting a hub gene's own primary perturbation over successive propagation iterations, since the delta matrix was fully overwritten by the network-propagated version at every step rather than persisting the hub signal ([#27](https://github.com/smorabit/compact/issues/27)). Hub gene rows are now reset to their initial perturbation value after every iteration (the "persistent-source hub model"), modeling a constitutive perturbation that holds hub genes at their new expression level while downstream genes respond. The signal-amplification warning is now correctly scoped to the non-hub subgraph.
- Fixed a related edge case in `PerturbationTransitions`: cells with an exact all-zero perturbation delta could get a spurious non-NaN transition weight to a specific neighbor (rather than staying in place), depending on internal behavior of the underlying correlation routine on degenerate inputs. Such cells' off-diagonal transition weights are now zeroed deterministically.

# 0.1.2 (08-06-2026)
## Changes
- Fixed a bug in `PredictPerturbationTime`, `PredictAttractors`, `PredictFates`, and `PredictCommitment` where cells with no valid outgoing transitions after graph masking produced an invalid sub-stochastic transition matrix — silently leaking probability mass, reporting false "already arrived" hitting times, and in some cases crashing `PredictAttractors`'s eigensolver on disconnected transition graphs ([#24](https://github.com/smorabit/compact/issues/24)). Zero-outdegree cells are now repaired with a self-transition instead of an invalid all-zero row.
- `PredictAttractors` now detects disconnected components in the transition graph before solving and restricts the eigensolve to the largest one, instead of potentially crashing.
- `FindConnectedComponents` gains a `matrix` argument to analyse a raw adjacency/transition matrix directly, in addition to a Seurat object graph.

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
