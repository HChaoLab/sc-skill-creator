# scanpy_leiden

![Python](https://img.shields.io/badge/language-Python-3776AB?logo=python)  ![Scanpy](https://img.shields.io/badge/package-Scanpy-lightgrey)

> Auto-generated from `skill.yaml`. Edit source files, not this README.

Leiden clustering partitions single cells into groups based on transcriptional similarity using community detection on a k-nearest neighbor graph. Identifies cell populations (cell types, states, subtypes) without prior labels. Resolution controls granularity: low resolution finds broad cell types, high resolution reveals subtypes and transitional states.

## Parameters

| Parameter | Default | Effect |
|-----------|---------|--------|
| `resolution` | `1.0` | Controls cluster granularity. Higher values produce more, smaller clusters. Too ... |
| `n_neighbors` | `15` | Determines k in kNN graph. More neighbors smooth transitions between populations... |
| `random_state` | `42` |  |
| `key_added` | `leiden` |  |

## Parameter Impact

### `resolution`

Controls cluster granularity. Higher values produce more, smaller clusters. Too few suggests over-aggregation; too many suggests noise is split as structure.

### `n_neighbors`

Determines k in kNN graph. More neighbors smooth transitions between populations. Fewer neighbors preserve local structure but may fragment clusters.

## Parameter Tuning Guide

The Agent uses this table to adjust parameters when critic feedback is received.

### `resolution`

| Condition | Adjust | Causal Chain |
|-----------|--------|--------------|
| `too_few_clusters` | increase by 0.3 | Higher resolution lowers modularity threshold, allowing more communiti... |
| `too_many_clusters` | decrease by 0.3 | Lower resolution raises modularity threshold, causing similar cells to... |
| `small_clusters` | decrease by 0.2 | Lower resolution merges small clusters with their nearest neighbors |

### `n_neighbors`

| Condition | Adjust | Causal Chain |
|-----------|--------|--------------|
| `noisy_fragmented` | increase by 7 | More neighbors smooth the kNN graph by averaging over more cells |
| `merged_undifferentiated` | decrease by 5 | Fewer neighbors sharpen neighborhood boundaries, revealing finer struc... |

## Critic Metrics

- `n_clusters`
- `min_cluster_size`
- `max_cluster_size`
- `mean_cluster_size`

## Context Parameters

Pass via `context` dict to `critic_post_process()`.

- `protocol`

## Error Handling

| Error | Solution |
|-------|----------|
| `MemoryError` | Reduce n_neighbors or n_pcs |
| `KeyError: neighbors` | Run sc.pp.neighbors() before sc.tl.leiden() |

## Execution Notes

scanpy_leiden — Leiden Clustering

PURPOSE:
    Partitions single cells into groups based on transcriptional similarity
    using community detection on a k-nearest neighbor graph.
    Resolution controls granularity: low values = broad cell types,
    high values = subtypes and transitional states.

INPUTS:
    input_data  : str | AnnData  — file path (.h5ad) or in-memory AnnData object
    params_dict : dict           — overrides defaults (optional)

DEFAULTS:
    resolution=1.0, n_neighbors=15, random_state=42, key_added='leiden'

OUTPUTS:
    AnnData with:
        adata.obs['leiden']           — cluster labels (string categories)
        adata.uns['analysis_history'] — execution record appended

## Critic Notes

scanpy_leiden critic — Evaluates Leiden clustering quality

Extracts hard metrics from the clustered AnnData and returns structured
feedback that maps directly to parameter_science_guide keys in skill.yaml.

FEEDBACK KEYS (must exactly match parameter_science_guide conditions):
    too_few_clusters    → n_clusters < threshold; resolution too low
    too_many_clusters   → n_clusters > threshold; resolution too high
    small_clusters      → min_cluster_size below threshold; marginal communities present
    noisy_fragmented    → high cluster count with many tiny clusters; n_neighbors too low
