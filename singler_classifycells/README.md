# singler_classifycells — SingleR Reference-Based Annotation

**Package**: SingleR + celldex | **Language**: R | **Function**: `SingleR::SingleR()`

## Purpose

Annotates cell types by correlating each cell's expression against curated reference datasets (celldex). Gold-standard automated method for well-characterized cell types. Accepts both Seurat and SingleCellExperiment objects.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.R` | `run_classifycells(input_data, params_list)` |
| `critic.R` | `critic_post_process(obj, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `ref` | `"HumanPrimaryCellAtlasData"` | celldex reference function name |
| `labels_field` | `"label.main"` | `label.main` = broad; `label.fine` = detailed |
| `de_method` | `"classic"` | Gene selection for scoring |
| `prune` | `TRUE` | Remove low-confidence assignments |

## Common Reference Datasets (celldex)

| Reference | Best For |
|-----------|----------|
| `HumanPrimaryCellAtlasData` | General human tissues |
| `BlueprintEncodeData` | Human immune & stromal |
| `ImmGenData` | Mouse immune cells |
| `MouseRNAseqData` | General mouse tissues |
| `DatabaseImmuneCellExpressionData` | Human PBMC/immune |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `too_few_cell_types` | Too coarse | Switch to `label.fine` |
| `too_many_pruned` | >20% cells pruned | Switch to `label.main` |
| `low_confidence_labels` | Low delta.next scores | Try a more specific reference |
| `high_na_rate` | Many NA labels | Set `prune = FALSE` |
