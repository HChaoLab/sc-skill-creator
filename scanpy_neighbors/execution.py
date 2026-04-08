"""
scanpy_neighbors — kNN Graph Construction

PURPOSE:
    Constructs a k-nearest neighbor graph in PCA (or corrected embedding) space.
    Foundation for Leiden/Louvain clustering and UMAP visualization.

INPUTS:
    input_data  : str | AnnData
    params_dict : dict (optional)

DEFAULTS:
    n_neighbors=15, n_pcs=40, use_rep='X_pca', metric='euclidean', random_state=42

OUTPUTS:
    AnnData with adata.uns['neighbors'], adata.obsp['connectivities'],
    adata.obsp['distances'], and adata.uns['analysis_history']
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


def run_neighbors(input_data, params_dict=None, default_params=None):
    # Memory-First
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)

    # Parameter Safety
    default_params = default_params or {
        'n_neighbors': 15,
        'n_pcs': 40,
        'use_rep': 'X_pca',
        'metric': 'euclidean',
        'random_state': 42,
    }
    current_params = {**default_params, **(params_dict or {})}

    # Validate use_rep exists
    use_rep = current_params['use_rep']
    if use_rep not in adata.obsm:
        raise KeyError(
            f"Representation '{use_rep}' not found in adata.obsm. "
            f"Available: {list(adata.obsm.keys())}. "
            "Run the corresponding embedding step first."
        )

    # Cap n_pcs to available dimensions
    n_pcs = min(current_params['n_pcs'], adata.obsm[use_rep].shape[1])

    sc.pp.neighbors(
        adata,
        n_neighbors=current_params['n_neighbors'],
        n_pcs=n_pcs,
        use_rep=use_rep,
        metric=current_params['metric'],
        random_state=current_params['random_state'],
    )

    # Collect metrics
    conn = adata.obsp['connectivities']
    mean_conn = float(conn.data.mean()) if conn.nnz > 0 else 0.0

    # Metadata Footprinting
    if 'analysis_history' not in adata.uns:
        adata.uns['analysis_history'] = []
    adata.uns['analysis_history'].append(_to_serializable({
        'skill_id': 'scanpy_neighbors',
        'params': {**current_params, 'n_pcs': n_pcs},
        'metrics': {
            'n_neighbors_used': current_params['n_neighbors'],
            'n_pcs_used': n_pcs,
            'use_rep': use_rep,
            'mean_connectivities': mean_conn,
            'n_cells': adata.n_obs,
        },
        'timestamp': datetime.now().isoformat()
    }))

    if not isinstance(input_data, ad.AnnData):
        adata.write(input_data)
    return adata
