"""
celltypist_annotate — Automated Cell Type Annotation

PURPOSE:
    Assigns probabilistic cell type labels using pre-trained CellTypist models.
    Majority voting refines predictions using cluster-level consensus.
    Requires log-normalized data (sc.pp.normalize_total + sc.pp.log1p).

INPUTS:
    input_data  : str | AnnData
    params_dict : dict (optional)

DEFAULTS:
    model='Immune_All_Low.pkl', majority_voting=True, over_clustering=None,
    min_prop=0.0, p_thres=0.5

OUTPUTS:
    AnnData with adata.obs['predicted_labels'], adata.obs['majority_voting']
    (if majority_voting=True), adata.obs['conf_score'],
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


def run_annotate(input_data, params_dict=None, default_params=None):
    import celltypist
    from celltypist import models

    # Memory-First
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)

    # Parameter Safety
    default_params = default_params or {
        'model': 'Immune_All_Low.pkl',
        'majority_voting': True,
        'over_clustering': None,
        'min_prop': 0.0,
        'p_thres': 0.5,
    }
    current_params = {**default_params, **(params_dict or {})}

    # Load model
    model_input = current_params['model']
    try:
        model = models.Model.load(model=model_input)
    except Exception:
        # Try downloading if not found locally
        models.download_models(force_update=False, model=model_input)
        model = models.Model.load(model=model_input)

    # Run annotation
    kwargs = {
        'adata': adata,
        'model': model,
        'majority_voting': current_params['majority_voting'],
        'p_thres': current_params['p_thres'],
        'min_prop': current_params['min_prop'],
    }
    over_clust = current_params.get('over_clustering')
    if over_clust is not None:
        kwargs['over_clustering'] = over_clust

    predictions = celltypist.annotate(**kwargs)

    # Write results back to adata
    adata = predictions.to_adata()

    # Collect metrics
    label_col = 'majority_voting' if current_params['majority_voting'] else 'predicted_labels'
    labels = adata.obs[label_col]
    n_cell_types = int(labels.nunique())
    pct_unassigned = float((labels == 'Unassigned').sum() / len(labels) * 100)

    conf_col = 'conf_score' if 'conf_score' in adata.obs.columns else None
    mean_max_prob = float(adata.obs[conf_col].mean()) if conf_col else float('nan')

    top3 = labels[labels != 'Unassigned'].value_counts().head(3).index.tolist()

    # Metadata Footprinting
    if 'analysis_history' not in adata.uns:
        adata.uns['analysis_history'] = []
    adata.uns['analysis_history'].append(_to_serializable({
        'skill_id': 'celltypist_annotate',
        'params': {
            'model': current_params['model'],
            'majority_voting': current_params['majority_voting'],
            'p_thres': current_params['p_thres'],
            'min_prop': current_params['min_prop'],
        },
        'metrics': {
            'n_cell_types': n_cell_types,
            'pct_unassigned': pct_unassigned,
            'mean_max_probability': mean_max_prob,
            'top3_celltypes': top3,
            'n_cells': adata.n_obs,
        },
        'timestamp': datetime.now().isoformat()
    }))

    if not isinstance(input_data, ad.AnnData):
        adata.write(input_data)
    return adata
