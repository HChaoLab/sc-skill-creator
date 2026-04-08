# sctype_score — ScType Marker Gene Scoring

**Package**: ScType (IanevskiAleksandr) | **Language**: R

## Purpose

Automatically assigns cell types to clusters by computing marker gene enrichment scores. Uses built-in tissue-specific databases or custom XLSX marker lists. Fully interpretable and customizable.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.R` | `run_sctype_score(input_data, params_list)` |
| `critic.R` | `critic_post_process(obj, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `tissue_type` | `"Immune system"` | See available tissues below |
| `db_path` | `NULL` | Path to custom XLSX; NULL = built-in DB |
| `scaled_assay` | `"scale.data"` | Use `"SCT"` if SCTransform was applied |
| `cluster_col` | `"seurat_clusters"` | Cluster column for aggregation |

## Available Built-in Tissues

`Immune system`, `Pancreas`, `Liver`, `Eye`, `Kidney`, `Brain`, `Lung`, `Adrenal`, `Heart`, `Intestine`, `Muscle`, `Placenta`, `Spleen`, `Stomach`, `Thymus`

## Custom Marker Database Format (XLSX)

| Column | Content |
|--------|---------|
| `tissueType` | Tissue name |
| `cellName` | Cell type name |
| `geneSymbolmore1` | Positive marker genes (comma-separated) |
| `geneSymbolmore2` | Negative marker genes (comma-separated) |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `high_unknown` | >30% Unknown clusters | Check `tissue_type`; provide custom DB |
| `low_scores` | Low enrichment scores | Verify `scaled_assay` slot is correct |
| `too_few_types` | Fewer types than expected | Try more specific tissue type |

## Prerequisites

`ScaleData()` or `SCTransform()` → clustering → **`sctype_score()`**
