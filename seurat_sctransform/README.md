# seurat_sctransform — SCTransform Normalization

**Package**: Seurat | **Language**: R | **Function**: `Seurat::SCTransform()`

## Purpose

Regularized negative binomial regression that removes sequencing depth confounding while preserving biological variance. Replaces the `NormalizeData → FindVariableFeatures → ScaleData` pipeline with a single, statistically principled step.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.R` | `run_sctransform(input_data, params_list)` |
| `critic.R` | `critic_post_process(obj, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `variable.features.n` | 3000 | HVF for PCA; range 1000–5000 |
| `vst.flavor` | `"v2"` | Recommended; use `"glmGamPoi"` for large datasets |
| `vars.to.regress` | `NULL` | e.g., `c("percent.mt", "S.Score", "G2M.Score")` |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `too_few_hvf` | Too few variable features | Increase `variable.features.n` |
| `noisy_features` | Too many variable features | Decrease `variable.features.n` |
| `cell_cycle_confounded` | Cell cycle driving clusters | Add `S.Score`/`G2M.Score` to `vars.to.regress` |
| `mito_confounded` | High mito content present | Add `percent.mt` to `vars.to.regress` |

## Prerequisites

Raw counts in `@assays$RNA@counts` → **`SCTransform()`** → `RunPCA()` → `FindNeighbors()`
