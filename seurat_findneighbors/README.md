# seurat_findneighbors — Seurat SNN Graph Construction

**Package**: Seurat | **Language**: R | **Function**: `Seurat::FindNeighbors()`

## Purpose

Constructs a Shared Nearest Neighbor (SNN) graph in PCA (or batch-corrected) embedding space. Direct prerequisite for `FindClusters()` and `RunUMAP()`.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.R` | `run_findneighbors(input_data, params_list)` |
| `critic.R` | `critic_post_process(obj, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `dims` | `1:20` | PCA dimensions to use |
| `k.param` | 20 | k for kNN; range 5–50 |
| `reduction` | `"pca"` | Use `"harmony"` after batch correction |
| `prune.SNN` | 0.0667 | Jaccard similarity pruning threshold |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `missing_rare_populations` | Too few dims | Increase `dims` range |
| `noisy_clusters` | Too many dims | Reduce `dims` range |
| `fragmented_clusters` | k.param too low | Increase `k.param` by 5 |
| `over_merged_populations` | k.param too high | Decrease `k.param` by 5 |

## Prerequisites

`RunPCA()` (or Harmony/scVI integration) → **`FindNeighbors()`** → `FindClusters()` + `RunUMAP()`
