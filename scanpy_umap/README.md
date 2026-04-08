# scanpy_umap — UMAP Embedding

**Package**: scanpy | **Language**: Python | **Function**: `sc.tl.umap`

## Purpose

Non-linear dimensionality reduction for 2D/3D visualization. Preserves local and global transcriptional structure. **Does not affect clustering** — used for visual interpretation only.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.py` | `run_umap(input_data, params_dict)` |
| `critic.py` | `critic_post_process(adata, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `min_dist` | 0.5 | Cluster compactness; range 0.0–1.0 |
| `spread` | 1.0 | Inter-cluster spacing scale |
| `n_components` | 2 | Use 3 for 3D visualization |
| `random_state` | 42 | Reproducibility seed |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `overlapping_clusters` | Clusters visually merged | Decrease `min_dist` |
| `too_compressed` | Embedding range too small | Increase `spread` |
| `clusters_too_close` | Low centroid separation | Increase `spread` |
| `clusters_too_spread` | Embedding over-dispersed | Decrease `spread` |

## Prerequisites

`sc.pp.neighbors()` → **`sc.tl.umap()`**
