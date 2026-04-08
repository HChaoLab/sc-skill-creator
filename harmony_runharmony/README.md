# harmony_runharmony — Harmony Batch Correction

**Package**: harmony | **Language**: R | **Function**: `harmony::RunHarmony()`

## Purpose

Corrects batch effects by iteratively aligning PCA embeddings across batches. Only modifies the corrected `harmony` reduction — original count matrices are untouched. After running, use `reduction = "harmony"` in `FindNeighbors()` and `RunUMAP()`.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.R` | `run_runharmony(input_data, params_list)` |
| `critic.R` | `critic_post_process(obj, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `group.by.vars` | `"orig.ident"` | Batch column(s) in metadata |
| `theta` | 2.0 | Diversity penalty; range 0–4 |
| `lambda` | 1.0 | Ridge penalty; range 0.5–2 |
| `sigma` | 0.1 | Soft-assignment width |
| `max.iter.harmony` | 10 | Convergence iterations |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `under_corrected` | Batches still separate | Increase `theta` by 0.5 |
| `over_corrected` | Biological structure lost | Decrease `theta` by 0.5; increase `lambda` |
| `convergence_not_reached` | Algorithm not converged | Increase `max.iter.harmony` by 5 |

## After Running

```r
# IMPORTANT: Use harmony reduction for all downstream steps
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)
```

## Prerequisites

`RunPCA()` (after `NormalizeData` or `SCTransform`) → **`RunHarmony()`** → `FindNeighbors(reduction="harmony")`
