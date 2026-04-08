"""
scanpy_leiden critic — Evaluates Leiden clustering quality

Extracts hard metrics from the clustered AnnData and returns structured
feedback that maps directly to parameter_science_guide keys in skill.yaml.

FEEDBACK KEYS (must exactly match parameter_science_guide conditions):
    too_few_clusters    → n_clusters < threshold; resolution too low
    too_many_clusters   → n_clusters > threshold; resolution too high
    small_clusters      → min_cluster_size below threshold; marginal communities present
    noisy_fragmented    → high cluster count with many tiny clusters; n_neighbors too low
"""


def critic_post_process(adata, context=None):
    """Evaluate clustering quality and return structured feedback.

    Parameters
    ----------
    adata : AnnData
        Clustered AnnData with adata.obs['leiden'] present.
    context : dict, optional
        Protocol context. Supported keys:
            protocol : str  — '10x Genomics' | 'Smart-seq2' | 'Drop-seq' | 'inDrop'

    Returns
    -------
    dict with keys:
        metrics      : dict  — extracted numeric measurements
        warnings     : list  — human-readable diagnostic messages
        feedback     : list  — condition strings matching parameter_science_guide keys
        success      : bool  — True when feedback is empty
        context_used : str   — which protocol thresholds were applied
    """
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')

    thresholds = {
        '10x Genomics': {'min_clusters': 3,  'max_clusters': 200, 'min_cluster_size': 5,  'noisy_ratio': 0.2},
        'Smart-seq2':   {'min_clusters': 5,  'max_clusters': 100, 'min_cluster_size': 10, 'noisy_ratio': 0.15},
        'Drop-seq':     {'min_clusters': 3,  'max_clusters': 150, 'min_cluster_size': 5,  'noisy_ratio': 0.2},
        'inDrop':       {'min_clusters': 3,  'max_clusters': 180, 'min_cluster_size': 5,  'noisy_ratio': 0.2},
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])

    metrics = {}
    warnings = []
    feedback = []

    # Guard: leiden must exist
    if 'leiden' not in adata.obs:
        return {
            'metrics': {},
            'warnings': ['CRITICAL: leiden column not found in adata.obs — was sc.tl.leiden() run?'],
            'feedback': [],
            'success': False,
            'context_used': protocol,
        }

    counts = adata.obs['leiden'].value_counts()
    n_clusters = adata.obs['leiden'].nunique()

    metrics['n_clusters']        = int(n_clusters)
    metrics['n_cells']           = int(adata.n_obs)
    metrics['min_cluster_size']  = int(counts.min())
    metrics['max_cluster_size']  = int(counts.max())
    metrics['mean_cluster_size'] = float(counts.mean())

    # ── Check 1: too few clusters ──
    if metrics['n_clusters'] < t['min_clusters']:
        warnings.append(
            f"UNDER-SEGMENTATION: only {metrics['n_clusters']} clusters "
            f"(expected >= {t['min_clusters']}) for {protocol}"
        )
        feedback.append('too_few_clusters')

    # ── Check 2: too many clusters ──
    if metrics['n_clusters'] > t['max_clusters']:
        warnings.append(
            f"OVER-SEGMENTATION: {metrics['n_clusters']} clusters "
            f"(expected <= {t['max_clusters']}) for {protocol}"
        )
        feedback.append('too_many_clusters')

    # ── Check 3: small clusters ──
    if metrics['min_cluster_size'] < t['min_cluster_size']:
        warnings.append(
            f"SMALL CLUSTERS: min cluster has {metrics['min_cluster_size']} cells "
            f"(expected >= {t['min_cluster_size']})"
        )
        feedback.append('small_clusters')

    # ── Check 4: noisy fragmentation ──
    # Many clusters where a high fraction are tiny → n_neighbors likely too low
    tiny = (counts < t['min_cluster_size']).sum()
    if metrics['n_clusters'] > 0:
        tiny_ratio = tiny / metrics['n_clusters']
        metrics['tiny_cluster_fraction'] = float(tiny_ratio)
        if tiny_ratio > t['noisy_ratio'] and 'small_clusters' not in feedback:
            warnings.append(
                f"NOISY FRAGMENTATION: {tiny} / {metrics['n_clusters']} clusters "
                f"are tiny ({tiny_ratio:.0%}); consider increasing n_neighbors"
            )
            feedback.append('noisy_fragmented')

    return {
        'metrics': metrics,
        'warnings': warnings,
        'feedback': feedback,    # Agent looks up each item in parameter_science_guide
        'success': len(feedback) == 0,
        'context_used': protocol,
    }
