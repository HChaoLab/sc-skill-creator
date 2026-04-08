"""
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
"""
import scanpy as sc
import anndata as ad
import numpy as np
from datetime import datetime


def _to_serializable(obj):
    """Recursively convert numpy types to Python natives for h5ad compatibility.

    adata.uns is saved to HDF5 when writing .h5ad. h5py rejects numpy scalars
    (np.int64, np.float32, etc.) and arrays — only Python native types are safe.
    Always wrap the analysis_history entry with this before appending.
    """
    if isinstance(obj, dict):
        return {str(k): _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(i) for i in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def run_leiden(input_data, params_dict=None, default_params=None):
    # Memory-First: accept both file path and in-memory AnnData
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)

    # Parameter Safety: merge caller params over defaults
    default_params = default_params or {
        'resolution': 1.0,
        'n_neighbors': 15,
        'random_state': 42,
        'key_added': 'leiden',
    }
    current_params = {**default_params, **(params_dict or {})}

    # Run Leiden clustering
    sc.tl.leiden(
        adata,
        resolution=current_params['resolution'],
        n_neighbors=current_params['n_neighbors'],
        random_state=current_params['random_state'],
        key_added=current_params['key_added'],
    )

    key = current_params['key_added']
    cluster_sizes = adata.obs[key].value_counts()

    # Metadata Footprinting — _to_serializable prevents h5ad TypeError on numpy types
    if 'analysis_history' not in adata.uns:
        adata.uns['analysis_history'] = []
    adata.uns['analysis_history'].append(_to_serializable({
        'skill_id': 'scanpy_leiden',
        'params': current_params,
        'metrics': {
            'n_clusters': adata.obs[key].nunique(),
            'n_cells': adata.n_obs,
            'min_cluster_size': cluster_sizes.min() if len(cluster_sizes) else 0,
            'max_cluster_size': cluster_sizes.max() if len(cluster_sizes) else 0,
        },
        'timestamp': datetime.now().isoformat(),
    }))

    # Write back only if input was a file path
    if not isinstance(input_data, ad.AnnData):
        adata.write(input_data)

    return adata
