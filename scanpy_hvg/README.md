# scanpy_hvg — Highly Variable Gene Selection

**Package**: scanpy | **Language**: Python | **Function**: `sc.pp.highly_variable_genes`

## Purpose

Identifies genes with greater cell-to-cell variation than expected by chance. These HVGs carry the most biological signal for distinguishing cell types and states, and are used to restrict downstream PCA, clustering, and embedding steps.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.py` | `run_hvg(input_data, params_dict)` |
| `critic.py` | `critic_post_process(adata, context)` |
| `README.md` | This file — auto-generated, do not edit manually |

## Default Parameters

| Parameter | Default | Range | Effect |
|-----------|---------|-------|--------|
| `n_top_genes` | 2000 | 500–5000 | Number of HVGs to select |
| `min_mean` | 0.0125 | 0.001–0.1 | Minimum mean expression |
| `max_mean` | 3.0 | 2–10 | Maximum mean expression (log scale) |
| `min_disp` | 0.5 | 0.1–2.0 | Minimum normalized dispersion |
| `flavor` | `'seurat'` | seurat / seurat_v3 / cell_ranger | HVG algorithm |
| `batch_key` | `None` | column name or null | Per-batch HVG selection |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `too_few_hvg` | Too few genes selected | Increase `n_top_genes` or decrease `min_disp` |
| `too_many_hvg` | Too many genes selected | Decrease `n_top_genes` |
| `low_variance_dominated` | Selected genes have low dispersion | Increase `min_disp` |
| `missing_rare_markers` | Very low pct of genes selected | Decrease `min_disp` |

## Prerequisites

- `sc.pp.normalize_total()` must be run before this step (for `flavor='seurat'` or `'cell_ranger'`)
- `sc.pp.log1p()` must be run before this step (for `flavor='seurat'` or `'cell_ranger'`)
- For `flavor='seurat_v3'`, use on raw counts (before normalization)
