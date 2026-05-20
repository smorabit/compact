# Package index

## All functions

- [`AlertHeatmap()`](https://smorabit.github.io/compact/reference/AlertHeatmap.md)
  : AlertHeatmap: visualize pathway-level perturbation metrics as a
  heatmap

- [`AlertSystemScore()`](https://smorabit.github.io/compact/reference/AlertSystemScore.md)
  : AlertSystem: score pathway gene sets for a given module perturbation

- [`ApplyPerturbation()`](https://smorabit.github.io/compact/reference/ApplyPerturbation.md)
  : Apply In-Silico Perturbation

- [`ApplyPropagation()`](https://smorabit.github.io/compact/reference/ApplyPropagation.md)
  : Propagate In-Silico Perturbation Signal Through a Gene Network

- [`BarplotShap()`](https://smorabit.github.io/compact/reference/BarplotShap.md)
  : Plot Bar Chart of Top SHAP Driver Genes

- [`BeeswarmplotShap()`](https://smorabit.github.io/compact/reference/BeeswarmplotShap.md)
  : Plot SHAP Beeswarm Summary for Top Driver Genes

- [`CheckSaturation()`](https://smorabit.github.io/compact/reference/CheckSaturation.md)
  : Check for Signal Saturation After Network Propagation

- [`CheckSignalDecay()`](https://smorabit.github.io/compact/reference/CheckSignalDecay.md)
  : Check for Propagation Signal Decay

- [`CombinePerturbAssays()`](https://smorabit.github.io/compact/reference/CombinePerturbAssays.md)
  : CombinePerturbAssays: create a composite perturb assay by delta
  superposition (sparse-safe; COUNTS-only)

- [`ComputeDistance()`](https://smorabit.github.io/compact/reference/ComputeDistance.md)
  : ComputeDistance

- [`CustomPerturbation()`](https://smorabit.github.io/compact/reference/CustomPerturbation.md)
  : CustomPerturbation

- [`FindConnectedComponents()`](https://smorabit.github.io/compact/reference/FindConnectedComponents.md)
  : FindConnectedComponents

- [`FindShapKeyDriver()`](https://smorabit.github.io/compact/reference/FindShapKeyDriver.md)
  : Identify Key Driver Genes Using SHAP and XGBoost

- [`HeatmapDistance()`](https://smorabit.github.io/compact/reference/HeatmapDistance.md)
  : HeatmapDistance: Generate Heatmaps for Original and Perturbed
  Matrices

- [`LoadAlertSystem()`](https://smorabit.github.io/compact/reference/LoadAlertSystem.md)
  : Load a previously saved AlertSystem state into a Seurat object

- [`LoadAlertSystemPathways()`](https://smorabit.github.io/compact/reference/LoadAlertSystemPathways.md)
  : Load AlertSystem pathway gene sets (local or GitHub)

- [`MacrostateTransitions()`](https://smorabit.github.io/compact/reference/MacrostateTransitions.md)
  : MacrostateTransitions

- [`ModelZINB()`](https://smorabit.github.io/compact/reference/ModelZINB.md)
  : ModelZINB

- [`ModulePerturbation()`](https://smorabit.github.io/compact/reference/ModulePerturbation.md)
  : ModulePerturbation

- [`PerturbationLog2FC()`](https://smorabit.github.io/compact/reference/PerturbationLog2FC.md)
  : Calculate Log2 Fold Change of In-Silico Perturbation

- [`PerturbationTransitions()`](https://smorabit.github.io/compact/reference/PerturbationTransitions.md)
  : PerturbationTransitions

- [`PerturbationVectors()`](https://smorabit.github.io/compact/reference/PerturbationVectors.md)
  : Calculate Perturbation Transition Vectors for Cell Embeddings

- [`PlotMacrostateTransitions()`](https://smorabit.github.io/compact/reference/PlotMacrostateTransitions.md)
  : PlotMacrostateTransitions

- [`PlotTransitionVectors()`](https://smorabit.github.io/compact/reference/PlotTransitionVectors.md)
  : Plot Transition Vectors on a Reduced Dimensional Embedding

- [`PredictAttractors()`](https://smorabit.github.io/compact/reference/PredictAttractors.md)
  : Predict Perturbation Attractor States

- [`PredictCommitment()`](https://smorabit.github.io/compact/reference/PredictCommitment.md)
  : Predict Fate Commitment (Committor Probabilities)

- [`PredictFates()`](https://smorabit.github.io/compact/reference/PredictFates.md)
  : Predict Perturbation Fates (Forward Diffusion)

- [`PredictPerturbationTime()`](https://smorabit.github.io/compact/reference/PredictPerturbationTime.md)
  : Predict Perturbation Pseudotime (Absorbing Markov Chain Hitting
  Time)

- [`ReloadShapOutput()`](https://smorabit.github.io/compact/reference/ReloadShapOutput.md)
  : Reload SHAP Output from Disk

- [`SampleZINB()`](https://smorabit.github.io/compact/reference/SampleZINB.md)
  : SampleZINB

- [`SaveAlertSystem()`](https://smorabit.github.io/compact/reference/SaveAlertSystem.md)
  : Save the full AlertSystem state to disk

- [`StampPerturbProvenance()`](https://smorabit.github.io/compact/reference/StampPerturbProvenance.md)
  : StampPerturbProvenance: record baseline provenance for a perturb
  assay

- [`TFPerturbation()`](https://smorabit.github.io/compact/reference/TFPerturbation.md)
  : TFPerturbation

- [`VectorFieldCoherence()`](https://smorabit.github.io/compact/reference/VectorFieldCoherence.md)
  : Compute Local Coherence of a Vector Field

- [`.create_distance_heatmap()`](https://smorabit.github.io/compact/reference/dot-create_distance_heatmap.md)
  : This function generates a heatmap from a given matrix. It handles
  both matrices and data frames (converting them to matrices if needed).
  It also allows customization of the color palette, axis visibility,
  and legend display.

- [`.edist()`](https://smorabit.github.io/compact/reference/dot-edist.md)
  : edist

- [`.eucldist()`](https://smorabit.github.io/compact/reference/dot-eucldist.md)
  : eucldist

- [`.get_upper_tri()`](https://smorabit.github.io/compact/reference/dot-get_upper_tri.md)
  :

  This function is used internally to extract the upper triangle of a
  matrix. It replaces the lower triangle of the input matrix with `NA`,
  which is useful when working with symmetrical matrices such as
  distance or correlation matrices where only the upper triangle is
  needed.

- [`.spearmandist()`](https://smorabit.github.io/compact/reference/dot-spearmandist.md)
  : spearmandist
