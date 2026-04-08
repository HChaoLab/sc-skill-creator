# seurat_labeltransfer — Seurat Label Transfer

**Package**: Seurat | **Language**: R | **Functions**: `FindTransferAnchors()` + `TransferData()`

## Purpose

Transfers cell type annotations from an annotated reference Seurat object to a query dataset using CCA anchor-based integration. Best used when a high-quality atlas is available for the same tissue.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.R` | `run_labeltransfer(input_data, params_list)` |
| `critic.R` | `critic_post_process(obj, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `reference_path` | **(required)** | Path to reference RDS Seurat object |
| `ref_label_col` | `"cell_type"` | Label column in reference@meta.data |
| `dims` | `1:30` | CCA dimensions for anchor finding |
| `normalization_method` | `"LogNormalize"` | Use `"SCT"` if both used SCTransform |
| `k_weight` | 100 | Neighbors for label weighting |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `low_prediction_score` | Mean score <0.5 | Increase `dims`; check reference tissue match |
| `noisy_transfer` | High score variance | Decrease `dims` |
| `fragmented_predictions` | Labels inconsistent in clusters | Increase `k_weight` |
| `over_smoothed` | Too few predicted types | Decrease `k_weight` |

## Usage Example

```r
result <- run_labeltransfer(query_obj, params_list = list(
    reference_path = "/path/to/reference.rds",
    ref_label_col  = "cell_type",
    dims           = 1:30
))
```
