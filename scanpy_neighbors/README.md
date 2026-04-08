# scanpy_neighbors — kNN Graph Construction

**Package**: scanpy | **Language**: Python | **Function**: `sc.pp.neighbors`

## Purpose

Constructs a k-nearest neighbor (kNN) graph in PCA (or batch-corrected) embedding space. This graph is the direct prerequisite for Leiden clustering and UMAP visualization.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.py` | `run_neighbors(input_data, params_dict)` |
| `critic.py` | `critic_post_process(adata, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `n_neighbors` | 15 | k for kNN; range 5–50 |
| `n_pcs` | 40 | PCA dims to use for distances |
| `use_rep` | `'X_pca'` | Use `'X_harmony'` after batch correction |
| `metric` | `'euclidean'` | Distance metric |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `fragmented_clusters` | n_neighbors too low | Increase `n_neighbors` by 5 |
| `over_merged_clusters` | n_neighbors too high | Decrease `n_neighbors` by 5 |
| `noisy_structure` | n_pcs too high | Decrease `n_pcs` by 10 |
| `missing_rare_populations` | n_pcs too low | Increase `n_pcs` by 10 |

## Prerequisites

`sc.tl.pca()` (or harmony/scVI integration) → **`sc.pp.neighbors()`** → `sc.tl.leiden()` + `sc.tl.umap()`
