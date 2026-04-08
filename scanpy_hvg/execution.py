"""
scanpy_hvg — Highly Variable Gene Selection

PURPOSE:
    Identifies genes with greater cell-to-cell variation than expected by chance.
    These genes carry biological signal for distinguishing cell types and states.
    Restricting downstream analysis to HVGs reduces noise and dimensionality.

INPUTS:
    input_data  : str | AnnData
    params_dict : dict (optional)

DEFAULTS:
    n_top_genes=2000, min_mean=0.0125, max_mean=3.0, min_disp=0.5,
    flavor='seurat', batch_key=None, subset=False

OUTPUTS:
    AnnData with adata.var['highly_variable'] and adata.uns['analysis_history']
"""
import scanpy as sc
import anndata as ad
import numpy as np
from datetime import datetime


def _to_serializable(obj):
    """Recursively convert numpy types to Python natives for h5ad compatibility.
    adata.uns is saved to HDF5; h5py rejects numpy scalars and arrays."""
    if isinstance(obj, dict):
        return {str(k): _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(i) for i in obj]
    if isinstance(obj, np.integer):   return int(obj)
    if isinstance(obj, np.floating):  return float(obj)
    if isinstance(obj, np.ndarray):   return obj.tolist()
    if isinstance(obj, np.bool_):     return bool(obj)
    return obj


def run_hvg(input_data, params_dict=None, default_params=None):
    # Memory-First
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)

    # Parameter Safety
    default_params = default_params or {
        'n_top_genes': 2000,
        'min_mean': 0.0125,
        'max_mean': 3.0,
        'min_disp': 0.5,
        'flavor': 'seurat',
        'batch_key': None,
        'subset': False,
    }
    current_params = {**default_params, **(params_dict or {})}

    # Build kwargs — exclude None values for optional args
    kwargs = {
        'n_top_genes': current_params['n_top_genes'],
        'min_mean': current_params['min_mean'],
        'max_mean': current_params['max_mean'],
        'min_disp': current_params['min_disp'],
        'flavor': current_params['flavor'],
        'subset': current_params['subset'],
    }
    if current_params.get('batch_key') is not None:
        kwargs['batch_key'] = current_params['batch_key']

    sc.pp.highly_variable_genes(adata, **kwargs)

    # Collect metrics
    n_hvg = int(adata.var['highly_variable'].sum())
    hvg_mask = adata.var['highly_variable']
    mean_disp_hvg = float(adata.var.loc[hvg_mask, 'dispersions_norm'].mean()) \
        if 'dispersions_norm' in adata.var.columns else float('nan')
    pct_hvg = float(n_hvg / adata.n_vars * 100)

    # Metadata Footprinting — _to_serializable prevents h5ad write errors
    if 'analysis_history' not in adata.uns:
        adata.uns['analysis_history'] = []
    adata.uns['analysis_history'].append(_to_serializable({
        'skill_id': 'scanpy_hvg',
        'params': current_params,
        'metrics': {
            'n_hvg': n_hvg,
            'mean_dispersion_hvg': mean_disp_hvg,
            'pct_hvg': pct_hvg,
            'n_genes_total': adata.n_vars,
        },
        'timestamp': datetime.now().isoformat()
    }))

    if not isinstance(input_data, ad.AnnData):
        adata.write(input_data)
    return adata
