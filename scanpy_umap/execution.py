"""
scanpy_umap — UMAP Embedding

PURPOSE:
    Non-linear dimensionality reduction for 2D/3D visualization of single-cell data.
    Preserves local and global transcriptional structure. Does NOT affect clustering.

INPUTS:
    input_data  : str | AnnData
    params_dict : dict (optional)

DEFAULTS:
    min_dist=0.5, spread=1.0, n_components=2, maxiter=None, random_state=42

OUTPUTS:
    AnnData with adata.obsm['X_umap'] and adata.uns['analysis_history']
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


def run_umap(input_data, params_dict=None, default_params=None):
    # Memory-First
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)

    # Parameter Safety
    default_params = default_params or {
        'min_dist': 0.5,
        'spread': 1.0,
        'n_components': 2,
        'maxiter': None,
        'random_state': 42,
    }
    current_params = {**default_params, **(params_dict or {})}

    kwargs = {
        'min_dist': current_params['min_dist'],
        'spread': current_params['spread'],
        'n_components': current_params['n_components'],
        'random_state': current_params['random_state'],
    }
    if current_params.get('maxiter') is not None:
        kwargs['maxiter'] = current_params['maxiter']

    sc.tl.umap(adata, **kwargs)

    # Collect metrics
    umap_coords = adata.obsm['X_umap']
    umap_range_x = float(umap_coords[:, 0].max() - umap_coords[:, 0].min())
    umap_range_y = float(umap_coords[:, 1].max() - umap_coords[:, 1].min())

    # Metadata Footprinting
    if 'analysis_history' not in adata.uns:
        adata.uns['analysis_history'] = []
    adata.uns['analysis_history'].append(_to_serializable({
        'skill_id': 'scanpy_umap',
        'params': current_params,
        'metrics': {
            'umap_range_x': umap_range_x,
            'umap_range_y': umap_range_y,
            'n_components': current_params['n_components'],
            'n_cells': adata.n_obs,
        },
        'timestamp': datetime.now().isoformat()
    }))

    if not isinstance(input_data, ad.AnnData):
        adata.write(input_data)
    return adata
