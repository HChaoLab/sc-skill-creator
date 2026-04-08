"""
scanpy_pca — Principal Component Analysis

PURPOSE:
    Linear dimensionality reduction that identifies orthogonal axes capturing
    maximum transcriptional variance. Foundational step before neighbor-graph
    construction, clustering, and UMAP visualization.

INPUTS:
    input_data  : str | AnnData
    params_dict : dict (optional)

DEFAULTS:
    n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state=42

OUTPUTS:
    AnnData with adata.obsm['X_pca'], adata.uns['pca'], adata.varm['PCs']
    and adata.uns['analysis_history']
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


def run_pca(input_data, params_dict=None, default_params=None):
    # Memory-First
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)

    # Parameter Safety
    default_params = default_params or {
        'n_comps': 50,
        'use_highly_variable': True,
        'svd_solver': 'arpack',
        'random_state': 42,
    }
    current_params = {**default_params, **(params_dict or {})}

    # If use_highly_variable requested but HVGs not computed, fall back gracefully
    use_hvg = current_params['use_highly_variable']
    if use_hvg and 'highly_variable' not in adata.var.columns:
        use_hvg = False

    sc.tl.pca(
        adata,
        n_comps=current_params['n_comps'],
        use_highly_variable=use_hvg,
        svd_solver=current_params['svd_solver'],
        random_state=current_params['random_state'],
    )

    # Collect metrics
    variance_ratio = adata.uns['pca']['variance_ratio']
    n_pcs_computed = len(variance_ratio)
    variance_ratio_top10 = float(np.sum(variance_ratio[:10]))
    cumulative_variance_50pcs = float(np.sum(variance_ratio[:min(50, n_pcs_computed)]))

    # Metadata Footprinting
    if 'analysis_history' not in adata.uns:
        adata.uns['analysis_history'] = []
    adata.uns['analysis_history'].append(_to_serializable({
        'skill_id': 'scanpy_pca',
        'params': {**current_params, 'use_highly_variable': use_hvg},
        'metrics': {
            'n_pcs_computed': n_pcs_computed,
            'variance_ratio_top10': variance_ratio_top10,
            'cumulative_variance_50pcs': cumulative_variance_50pcs,
            'variance_ratio_pc1': float(variance_ratio[0]),
            'variance_ratio_pc2': float(variance_ratio[1]) if n_pcs_computed > 1 else 0.0,
        },
        'timestamp': datetime.now().isoformat()
    }))

    if not isinstance(input_data, ad.AnnData):
        adata.write(input_data)
    return adata
