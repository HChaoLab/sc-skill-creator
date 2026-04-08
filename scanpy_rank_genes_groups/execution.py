"""
scanpy_rank_genes_groups — Differential Expression / Marker Gene Detection

PURPOSE:
    Identifies genes that distinguish each cell cluster from all others (one-vs-rest).
    Essential for cell type annotation and biological interpretation of clusters.

INPUTS:
    input_data  : str | AnnData
    params_dict : dict (optional)

DEFAULTS:
    groupby='leiden', method='wilcoxon', n_genes=100, use_raw=None,
    pts=True, key_added='rank_genes_groups'

OUTPUTS:
    AnnData with adata.uns['rank_genes_groups'] and adata.uns['analysis_history']
"""
import scanpy as sc
import anndata as ad
import numpy as np
from datetime import datetime


def _to_serializable(obj):
    """Recursively convert numpy types to Python natives for h5ad compatibility."""
    if isinstance(obj, dict):
        return {str(k): _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(i) for i in obj]
    if isinstance(obj, np.integer):   return int(obj)
    if isinstance(obj, np.floating):  return float(obj)
    if isinstance(obj, np.ndarray):   return obj.tolist()
    if isinstance(obj, np.bool_):     return bool(obj)
    return obj


def run_rank_genes_groups(input_data, params_dict=None, default_params=None):
    # Memory-First
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)

    # Parameter Safety
    default_params = default_params or {
        'groupby': 'leiden',
        'method': 'wilcoxon',
        'n_genes': 100,
        'use_raw': None,
        'pts': True,
        'key_added': 'rank_genes_groups',
    }
    current_params = {**default_params, **(params_dict or {})}

    # Auto-detect groupby fallback
    groupby = current_params['groupby']
    if groupby not in adata.obs.columns:
        # Try common alternatives
        for alt in ['louvain', 'leiden', 'cluster', 'cell_type']:
            if alt in adata.obs.columns:
                groupby = alt
                break
        else:
            raise KeyError(
                f"Column '{current_params['groupby']}' not found in adata.obs. "
                f"Available: {list(adata.obs.columns)}. Run clustering first."
            )

    # Auto-detect use_raw
    use_raw = current_params['use_raw']
    if use_raw is None:
        use_raw = adata.raw is not None

    kwargs = {
        'groupby': groupby,
        'method': current_params['method'],
        'n_genes': current_params['n_genes'],
        'use_raw': use_raw,
        'pts': current_params['pts'],
        'key_added': current_params['key_added'],
    }

    sc.tl.rank_genes_groups(adata, **kwargs)

    # Collect metrics
    key = current_params['key_added']
    result = adata.uns[key]
    groups = list(result['names'].dtype.names)
    n_groups = len(groups)

    # Top-1 score per group
    top1_scores = []
    for g in groups:
        scores = result['scores'][g]
        if len(scores) > 0:
            top1_scores.append(float(scores[0]))

    mean_top1_score = float(np.mean(top1_scores)) if top1_scores else 0.0
    min_top1_score = float(np.min(top1_scores)) if top1_scores else 0.0
    mean_genes_per_group = float(np.mean([
        len(result['names'][g]) for g in groups
    ]))

    # Metadata Footprinting
    if 'analysis_history' not in adata.uns:
        adata.uns['analysis_history'] = []
    adata.uns['analysis_history'].append(_to_serializable({
        'skill_id': 'scanpy_rank_genes_groups',
        'params': {**current_params, 'groupby': groupby, 'use_raw': use_raw},
        'metrics': {
            'n_groups': n_groups,
            'mean_top1_score': mean_top1_score,
            'min_top1_score': min_top1_score,
            'mean_genes_per_group': mean_genes_per_group,
        },
        'timestamp': datetime.now().isoformat()
    }))

    if not isinstance(input_data, ad.AnnData):
        adata.write(input_data)
    return adata
