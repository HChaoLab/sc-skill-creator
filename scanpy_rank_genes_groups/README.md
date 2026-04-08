# scanpy_rank_genes_groups — Marker Gene Detection

**Package**: scanpy | **Language**: Python | **Function**: `sc.tl.rank_genes_groups`

## Purpose

Identifies genes that distinguish each cell cluster from all others (one-vs-rest differential expression). Essential for cell type annotation and biological interpretation.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.py` | `run_rank_genes_groups(input_data, params_dict)` |
| `critic.py` | `critic_post_process(adata, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `groupby` | `'leiden'` | Cluster column in `adata.obs` |
| `method` | `'wilcoxon'` | Most robust for scRNA-seq |
| `n_genes` | 100 | Top genes per group |
| `pts` | `True` | Compute expression fractions |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `too_few_markers` | Not enough marker genes returned | Increase `n_genes` |
| `low_specificity` | Marker scores are weak | Check clustering; try `method='wilcoxon'` |
| `high_pct_background` | Markers expressed in most clusters | Re-cluster or adjust thresholds |

## Prerequisites

`sc.tl.leiden()` or any clustering → **`sc.tl.rank_genes_groups()`**
