"""
scanpy_pca critic — Evaluates PCA quality

FEEDBACK KEYS (must match parameter_science_guide conditions in skill.yaml):
    insufficient_variance_captured → top-10 PCs explain too little variance
    noise_dominated                → PC1 variance ratio very high (data may be dominated by one source)
"""


def critic_post_process(adata, context=None):
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')

    thresholds = {
        '10x Genomics': {'min_var_ratio_top10': 0.30, 'max_var_ratio_pc1': 0.60},
        'Smart-seq2':   {'min_var_ratio_top10': 0.40, 'max_var_ratio_pc1': 0.70},
        'Drop-seq':     {'min_var_ratio_top10': 0.25, 'max_var_ratio_pc1': 0.60},
        'CITE-seq':     {'min_var_ratio_top10': 0.25, 'max_var_ratio_pc1': 0.55},
        'Multiome':     {'min_var_ratio_top10': 0.20, 'max_var_ratio_pc1': 0.55},
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])

    metrics = {}
    warnings = []
    feedback = []

    if 'X_pca' not in adata.obsm:
        return {
            'metrics': {},
            'warnings': ['CRITICAL: X_pca not in adata.obsm — PCA step may not have run'],
            'feedback': [],
            'success': False,
            'context_used': protocol
        }

    import numpy as np
    variance_ratio = adata.uns['pca']['variance_ratio']
    n_pcs = len(variance_ratio)

    metrics['n_pcs_computed'] = n_pcs
    metrics['variance_ratio_top10'] = float(np.sum(variance_ratio[:min(10, n_pcs)]))
    metrics['variance_ratio_pc1'] = float(variance_ratio[0])
    metrics['cumulative_variance_50pcs'] = float(np.sum(variance_ratio[:min(50, n_pcs)]))

    if metrics['variance_ratio_top10'] < t['min_var_ratio_top10']:
        warnings.append(
            f"LOW VARIANCE: top-10 PCs explain only {metrics['variance_ratio_top10']:.1%} "
            f"(threshold {t['min_var_ratio_top10']:.1%} for {protocol}). "
            "Consider increasing n_comps or checking normalization."
        )
        feedback.append('insufficient_variance_captured')

    if metrics['variance_ratio_pc1'] > t['max_var_ratio_pc1']:
        warnings.append(
            f"DOMINATED BY PC1: PC1 explains {metrics['variance_ratio_pc1']:.1%} of variance "
            f"(threshold {t['max_var_ratio_pc1']:.1%}). "
            "Possible batch effect or sequencing depth dominance — check regress_out."
        )
        feedback.append('noise_dominated')

    return {
        'metrics': metrics,
        'warnings': warnings,
        'feedback': feedback,
        'success': not bool(feedback),
        'context_used': protocol
    }
