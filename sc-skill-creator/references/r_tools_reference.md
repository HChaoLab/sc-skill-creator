# R Single-Cell Analysis Tools — Biological Reference

This reference provides biological context for **commonly used** R-based single-cell analysis packages. It is not exhaustive — sc-skill-creator can generate skills for any R package. Use these entries as a guide when writing `cognitive_layer.purpose`, `parameter_impact`, and `parameter_science_guide` for any tool.

---

## Seurat

**Object**: `Seurat` (S4) — stores raw counts in `@assays$RNA@counts`, normalized data in `@assays$RNA@data`, embeddings in `@reductions`, cluster labels in `@meta.data`.

### SCTransform (normalization + HVG)
**Biological role**: Regularized negative binomial regression to remove sequencing depth confounding while preserving biological variance. Replaces the `NormalizeData → ScaleData → FindVariableFeatures` pipeline.

Key parameters:
- `variable.features.n` (default 3000): Number of variable features for PCA. More genes = more biological signal but slower.
- `vars.to.regress`: Covariates to regress out (e.g., `percent.mt`, `S.Score`, `G2M.Score`). Regressing cell cycle removes batch-like variation but may mask genuine cycling populations.
- `vst.flavor`: `"v2"` is recommended; uses improved Pearson residuals.

Critic metrics: variance explained by top PCs, count-depth correlation after normalization, % mitochondrial correlation.

### FindNeighbors (SNN graph construction)
**Biological role**: Constructs a shared nearest neighbor (SNN) graph in PCA space. Prerequisite for `FindClusters` and `RunUMAP`.

Key parameters:
- `dims` (default `1:10`): PCA dimensions used. Too few = miss structure; too many = noise dominates.
- `k.param` (default 20): k in kNN. Larger k smooths transitions, smaller k preserves local structure.
- `reduction` (default `"pca"`): Can use `"harmony"` or `"scvi"` after integration.

### FindClusters (community detection)
**Biological role**: Applies Louvain or Leiden algorithm on the SNN graph to identify cell populations. High resolution → more, smaller clusters (subtypes); low resolution → fewer, broader clusters (major cell types).

Key parameters:
- `resolution` (default 0.8): Primary granularity knob.
- `algorithm`: 1=Louvain, 2=Louvain with multilevel refinement, 3=SLM, 4=Leiden (more stable, recommended).
- `random.seed`: For reproducibility.

Feedback conditions: `too_few_clusters`, `too_many_clusters`, `small_clusters`, `resolution_plateau`.

### RunPCA
**Biological role**: Linear dimensionality reduction capturing maximum transcriptional variance axes. First step before neighbor graph and UMAP.

Key parameters:
- `npcs` (default 50): Number of PCs to compute. Use ElbowPlot to select meaningful dims.
- `features`: Subset to HVGs or custom gene list.

### RunUMAP
**Biological role**: Non-linear dimensionality reduction for 2D/3D visualization. Does NOT affect clustering; used for interpretation only.

Key parameters:
- `dims`: Same PCs used in `FindNeighbors` (must match for consistency).
- `n.neighbors` (default 30): Controls local vs global structure balance.
- `min.dist` (default 0.3): Controls cluster compactness. Lower = tighter clusters.
- `spread` (default 1): Controls inter-cluster spacing.
- `metric` (default `"cosine"`): Distance metric in high-dimensional space.

### FindMarkers / FindAllMarkers (differential expression)
**Biological role**: Identifies genes that distinguish one cell population from others. Used for cell type annotation.

Key parameters:
- `test.use`: Statistical test — `"wilcox"` (default, robust), `"bimod"` (bimodal), `"DESeq2"` (pseudo-bulk), `"MAST"` (hurdle model).
- `min.pct` (default 0.1): Minimum fraction of cells expressing gene. Filters lowly expressed genes.
- `logfc.threshold` (default 0.25): Minimum log fold-change. Reduces false positives.
- `only.pos` (default FALSE): Return only upregulated markers.

Feedback conditions: `too_few_markers`, `low_specificity`, `high_pct_background`.

### IntegrateLayers / IntegrateData (batch correction)
**Biological role**: Corrects batch effects while preserving biological variation. Multiple methods available.

Key parameters:
- `method`: `HarmonyIntegration`, `CCAIntegration`, `RPCAIntegration`, `scVIIntegration`.
- `normalization.method`: `"LogNormalize"` or `"SCT"`.
- `dims.to.integrate` / `k.weight`: Controls integration aggressiveness.

---

## Harmony (harmony / SeuratWrappers)

**Object**: Operates on PCA embedding matrix. Returns corrected PCA embedding; does NOT modify expression values.

### RunHarmony
**Biological role**: Iterative clustering in PCA space to align batch-specific distributions. Fast and memory-efficient. Works on `Seurat`, `SingleCellExperiment`, or raw matrix.

Key parameters:
- `group.by.vars` (required): Column(s) in metadata that encode batch (e.g., `"orig.ident"`, `c("donor", "technology")`).
- `theta` (default 2): Diversity penalty — higher values enforce more mixing across batches. Lower values allow biological substructure to partially resist merging.
- `lambda` (default 1): Ridge regression penalty on cluster centroids. Higher = more conservative correction.
- `sigma` (default 0.1): Width of soft-assignment kernel. Larger = softer assignments.
- `nclust` (default NULL, auto): Number of Harmony clusters. Auto-selected based on cell count.
- `max.iter.harmony` (default 10): Convergence iterations.

Feedback conditions: `over_corrected` (biological structure lost), `under_corrected` (batches still separate), `convergence_not_reached`.

**After RunHarmony**: Use `reduction = "harmony"` in `FindNeighbors` and `RunUMAP`.

---

## DESeq2

**Object**: `DESeqDataSet` (S4). Requires count matrix (integer), colData (sample metadata), and design formula.

### DESeq (main pipeline)
**Biological role**: Estimates size factors, dispersions, and fits negative binomial GLM. Gold standard for bulk and pseudo-bulk differential expression.

Key parameters:
- `test` (default `"Wald"`): `"Wald"` for pairwise comparison; `"LRT"` for multi-factor / interaction models.
- `fitType` (default `"parametric"`): Dispersion estimation method. `"local"` for datasets with few samples; `"mean"` for very few replicates.
- `useT` (default FALSE): Use t-distribution for small sample sizes (<6 per group).
- `minReplicatesForReplace` (default 7): Samples needed before replacing outlier counts.

Critic metrics: `% genes with padj < 0.05`, dispersion-mean relationship quality, size factor range.

Feedback conditions: `too_few_de_genes`, `dispersion_fit_poor`, `size_factor_extreme`.

### results / lfcShrink (effect size)
**Biological role**: `lfcShrink` applies adaptive shrinkage to log fold-changes, reducing noise for lowly expressed genes. Essential before ranking or visualization.

Key parameters:
- `type`: `"apeglm"` (recommended, adaptive), `"ashr"` (more flexible), `"normal"` (legacy).
- `alpha` (default 0.1): FDR threshold for `summary()`.

---

## edgeR

**Object**: `DGEList`. Stores count matrix, sample metadata, normalization factors.

### glmQLFit + glmQLFTest (quasi-likelihood pipeline)
**Biological role**: Quasi-likelihood F-test accounts for gene-specific variance inflation. Recommended for scRNA-seq pseudo-bulk and well-designed bulk experiments.

Key parameters:
- `robust` (default FALSE in older versions): Robust dispersion estimation — set TRUE for outlier samples.
- `abundance.trend` (default TRUE): Account for mean-abundance trend in dispersions.

### calcNormFactors
**Biological role**: Computes TMM (trimmed mean of M-values) normalization factors to correct for library composition bias.

Key parameters:
- `method` (default `"TMM"`): `"TMM"`, `"TMMwsp"` (with singleton pairing, better for sparse data), `"RLE"`, `"upperquartile"`, `"none"`.

---

## MAST

**Object**: `SingleCellAssay` or accepts Seurat/SCE objects via conversion.

### zlm (hurdle model)
**Biological role**: Two-part hurdle model — logistic component for dropout (zero/non-zero) and Gaussian component for expression level. Appropriate for scRNA-seq where many zeros are biological.

Key parameters:
- `formula`: Model formula, e.g., `~condition + cngeneson` (cellular detection rate as covariate).
- `method` (default `"bayesglm"`): `"bayesglm"` (Bayesian GLM, regularized), `"glm"`, `"glmer"`.
- `ebayes` (default TRUE): Empirical Bayes variance shrinkage.
- `ebayesControl`: List with `method` and `model.matrix` parameters for shrinkage.

Feedback conditions: `too_few_hurdle_de`, `model_convergence_issues`, `cngeneson_not_included`.

---

## Monocle3

**Object**: `cell_data_set` (CDS). Stores expression in sparse matrix, cell metadata in `colData`, gene metadata in `rowData`, UMAP in `reducedDims`.

### learn_graph (trajectory graph)
**Biological role**: Fits a principal graph through the UMAP embedding, representing developmental or activation trajectories.

Key parameters:
- `use_partition` (default TRUE): Learn separate graphs for each UMAP partition. Set FALSE to connect partitions into one trajectory.
- `close_loop` (default TRUE): Allow circular trajectories (e.g., cell cycle).
- `learn_graph_control`: List with `minimal_branch_len`, `geodesic_distance_ratio`, etc.

### order_cells (pseudotime)
**Biological role**: Assigns pseudotime values (developmental position) to each cell by projecting onto the principal graph from a user-specified root.

Key parameters:
- `root_cells`: Manually specify root cell IDs (e.g., progenitor cells). If NULL, interactive selection.
- `root_pr_nodes`: Alternatively specify root principal graph node.

Feedback conditions: `trajectory_too_branched`, `disconnected_trajectory`, `pseudotime_not_monotone`.

### cluster_cells (Leiden clustering for Monocle3)
**Biological role**: Groups cells for trajectory partitioning. Uses Leiden algorithm on UMAP neighbors.

Key parameters:
- `resolution` (default NULL, auto): Cluster resolution. Set explicitly for reproducibility.
- `k` (default 20): kNN for graph construction.
- `partition_qval` (default 0.05): FDR threshold for identifying separate trajectory partitions.

---

## scater / scran

**Object**: `SingleCellExperiment` (SCE). Stores counts in `assays$counts`, normalized in `assays$logcounts`, metadata in `colData`/`rowData`.

### perCellQCMetrics / addPerCellQCMetrics (QC)
**Biological role**: Computes per-cell QC metrics: total counts, detected genes, % mitochondrial, % ERCC spike-in.

Key parameters:
- `subsets`: Named list of gene subsets (e.g., `list(Mito = mito_genes)`).
- `use_altexps`: Include alternative experiments (ERCC, CITE-seq).

### quickPerCellQC (filtering)
**Biological role**: Removes outlier cells based on 3 MAD (median absolute deviation) rule applied to QC metrics. Protocol-aware thresholds.

Key parameters:
- `sum.field`, `detected.field`, `subsets_Mito_percent.field`: Which metrics to filter on.
- `nmads` (default 3): Number of MADs for outlier definition. Reduce to 2.5 for stricter QC.
- `type` (default `"lower"`): `"lower"` for count/gene thresholds; `"higher"` for mito %.

### computeSumFactors (scran normalization)
**Biological role**: Deconvolution normalization — pools cells to estimate size factors, then deconvolves. Handles near-zero counts better than simple library size normalization.

Key parameters:
- `clusters` (optional): Pre-computed cluster assignments to pool within clusters. Recommended for heterogeneous datasets.
- `min.mean` (default 0.1): Minimum mean count for genes used in pooling.

### modelGeneVar / getTopHVGs (HVG selection)
**Biological role**: Decomposes per-gene variance into technical (Poisson) and biological components. `getTopHVGs` returns top biologically variable genes.

Key parameters:
- `block`: Blocking factor for multi-sample datasets (models variance per sample).
- `span` (default 0.3): Loess smoothing span for the mean-variance trend.
- `n` in `getTopHVGs`: Number of HVGs to return (default 2000).

---

## Common R Robustness Patterns

### Memory-First (all R object types)
```r
if (is.character(input_data)) {
    obj <- readRDS(input_data)
} else {
    obj <- input_data
}
```

### Parameter Safety
```r
default_params <- list(resolution = 0.8, k.param = 20L)
current_params <- modifyList(default_params,
                             if (!is.null(params_list)) params_list else list())
```

### Metadata Footprinting — Seurat
```r
if (is.null(obj@misc$analysis_history)) obj@misc$analysis_history <- list()
obj@misc$analysis_history <- c(obj@misc$analysis_history, list(list(
    skill_id  = "seurat_findclusters",
    params    = current_params,
    metrics   = list(n_clusters = n_clusters),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)))
```

### Metadata Footprinting — SingleCellExperiment
```r
if (is.null(S4Vectors::metadata(sce)$analysis_history))
    S4Vectors::metadata(sce)$analysis_history <- list()
S4Vectors::metadata(sce)$analysis_history <- c(
    S4Vectors::metadata(sce)$analysis_history,
    list(list(
        skill_id  = "scran_cluster",
        params    = current_params,
        metrics   = list(n_clusters = n_clusters),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    ))
)
```

### Context-Aware Thresholds
```r
critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol <- if (!is.null(context$protocol)) context$protocol else "10x Genomics"

    thresholds <- list(
        "10x Genomics" = list(min_clusters = 3L, max_clusters = 200L),
        "Smart-seq2"   = list(min_clusters = 5L, max_clusters = 100L)
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]
    # use t$min_clusters, t$max_clusters, etc.
}
```

---

## Protocol Context Reference

| Protocol | Typical cells | Typical genes detected | Notes |
|----------|--------------|------------------------|-------|
| 10x Genomics (3') | 500–50,000 | 1,000–4,000 | High dropout, UMI-based |
| 10x Genomics (5') | 500–20,000 | 1,500–5,000 | Better for T/B cell receptors |
| Smart-seq2 | 100–10,000 | 4,000–10,000 | Full-length, plate-based, low dropout |
| Drop-seq | 1,000–20,000 | 800–3,000 | Similar to 10x but older |
| inDrop | 500–15,000 | 1,000–3,000 | Similar to 10x |
| CITE-seq | 1,000–20,000 | 1,000–4,000 | RNA + protein |
| Multiome (10x) | 1,000–20,000 | 2,000–6,000 | RNA + ATAC |

Use these as guide values when setting `context_params` and default thresholds.
