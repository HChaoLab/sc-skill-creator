"""
scanpy_neighbors critic — Evaluates kNN graph quality

FEEDBACK KEYS (must match parameter_science_guide conditions in skill.yaml):
    fragmented_clusters    → n_neighbors too low; graph under-connected
    over_merged_clusters   → n_neighbors too high; neighborhoods too large
    noisy_structure        → n_pcs too high; noise in distance calculation
    missing_rare_populations → n_pcs too low; rare variation excluded
"""
import numpy as np


def critic_post_process(adata, context=None):
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')

    thresholds = {
        '10x Genomics': {'min_neighbors': 5,  'max_neighbors': 50, 'min_conn': 0.01, 'max_conn': 0.5},
        'Smart-seq2':   {'min_neighbors': 10, 'max_neighbors': 50, 'min_conn': 0.02, 'max_conn': 0.6},
        'Drop-seq':     {'min_neighbors': 5,  'max_neighbors': 40, 'min_conn': 0.01, 'max_conn': 0.5},
        'CITE-seq':     {'min_neighbors': 10, 'max_neighbors': 50, 'min_conn': 0.01, 'max_conn': 0.5},
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])

    metrics = {}
    warnings = []
    feedback = []

    if 'neighbors' not in adata.uns or 'connectivities' not in adata.obsp:
        return {
            'metrics': {},
            'warnings': ['CRITICAL: neighbors not computed — run sc.pp.neighbors() first'],
            'feedback': [],
            'success': False,
            'context_used': protocol
        }

    params_used = adata.uns['neighbors'].get('params', {})
    n_neighbors = params_used.get('n_neighbors', 0)
    conn = adata.obsp['connectivities']
    mean_conn = float(conn.data.mean()) if conn.nnz > 0 else 0.0

    metrics['n_neighbors_used'] = n_neighbors
    metrics['mean_connectivities'] = mean_conn
    metrics['graph_density'] = float(conn.nnz / (adata.n_obs ** 2))

    if n_neighbors < t['min_neighbors']:
        warnings.append(
            f"UNDER-CONNECTED GRAPH: n_neighbors={n_neighbors} is very low "
            f"(min {t['min_neighbors']} for {protocol}). Graph may be fragmented."
        )
        feedback.append('fragmented_clusters')

    if n_neighbors > t['max_neighbors']:
        warnings.append(
            f"OVER-CONNECTED GRAPH: n_neighbors={n_neighbors} is high "
            f"(max {t['max_neighbors']} for {protocol}). May over-merge distinct populations."
        )
        feedback.append('over_merged_clusters')

    if mean_conn < t['min_conn']:
        warnings.append(
            f"WEAK CONNECTIVITY: mean connectivity={mean_conn:.4f} is very low. "
            "Graph may not reflect biological similarity well."
        )
        feedback.append('noisy_structure')

    return {
        'metrics': metrics,
        'warnings': warnings,
        'feedback': feedback,
        'success': not bool(feedback),
        'context_used': protocol
    }
