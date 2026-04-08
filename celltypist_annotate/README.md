# celltypist_annotate — Automated Cell Type Annotation

**Package**: celltypist | **Language**: Python | **Function**: `celltypist.annotate()`

## Purpose

Assigns probabilistic cell type labels to each cell using logistic regression models trained on curated single-cell atlases. Majority voting refines predictions using cluster-level consensus. Best suited for immune cell annotation.

## Files

| File | Description |
|------|-------------|
| `skill.yaml` | Metadata, cognitive layer, parameter science guide |
| `execution.py` | `run_annotate(input_data, params_dict)` |
| `critic.py` | `critic_post_process(adata, context)` |
| `README.md` | This file |

## Default Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `model` | `'Immune_All_Low.pkl'` | Pre-trained model; use `celltypist.models.models_description()` to list all |
| `majority_voting` | `True` | Refine with cluster consensus |
| `p_thres` | `0.5` | Probability threshold; lower = more labels |
| `min_prop` | `0.0` | Min fraction for majority voting |

## Common Models

| Model | Use Case |
|-------|----------|
| `Immune_All_Low.pkl` | Major immune cell types |
| `Immune_All_High.pkl` | Immune cell subtypes |
| `Healthy_COVID19_PBMC.pkl` | PBMC with COVID-19 context |
| `Human_Lung_Atlas.pkl` | Lung cell types |

## Critic Feedback Keys

| Key | Meaning | Suggested Fix |
|-----|---------|---------------|
| `high_unassigned` | >20% cells unassigned | Decrease `p_thres`; try different model |
| `low_confidence_labels` | Mean probability <0.55 | Model may not match tissue; try tissue-specific model |
| `over_homogenized` | Too few cell types detected | Decrease `min_prop` |
| `fragmented_voting` | Many mixed clusters unassigned | Decrease `min_prop` |

## Prerequisites

`sc.pp.normalize_total()` → `sc.pp.log1p()` → clustering (leiden) → **`celltypist.annotate()`**
