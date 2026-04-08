# scanpy_pca — Principal Component Analysis

**Package**: scanpy | **Language**: Python | **Function**: `sc.tl.pca`

## Purpose

Linear dimensionality reduction that identifies orthogonal axes (PCs) capturing maximum transcriptional variance. Reduces the high-dimensional gene space to 10–50 PCs before neighbor-graph, clustering, and UMAP steps.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.py` | `run_pca(input_data, params_dict)` |
| `critic.py` | `critic_post_process(adata, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `n_comps` | 50 | Number of PCs to compute |
| `use_highly_variable` | `True` | Restrict to HVGs if available |
| `svd_solver` | `'arpack'` | Use `'randomized'` for large datasets |
| `random_state` | 42 | Reproducibility seed |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `insufficient_variance_captured` | Top-10 PCs explain <30% variance | Increase `n_comps`; check normalization |
| `noise_dominated` | PC1 explains >60% variance | Run `sc.pp.regress_out` for depth/MT |

## Prerequisites

- `sc.pp.normalize_total()` → `sc.pp.log1p()` → `sc.pp.highly_variable_genes()` → **`sc.tl.pca()`**
