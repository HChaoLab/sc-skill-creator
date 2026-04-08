"""
scanpy_hvg critic — Evaluates HVG selection quality

FEEDBACK KEYS (must match parameter_science_guide conditions in skill.yaml):
    too_few_hvg          → n_top_genes or min_disp too restrictive
    too_many_hvg         → n_top_genes too high; noise included
    low_variance_dominated → mean dispersion of selected HVGs is low
    missing_rare_markers → pct_hvg very low relative to dataset size
"""


def critic_post_process(adata, context=None):
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')
    n_cells = context.get('n_cells', adata.n_obs)

    thresholds = {
        '10x Genomics': {'min_hvg': 500,  'max_hvg': 5000, 'min_disp_mean': 0.4},
        'Smart-seq2':   {'min_hvg': 1000, 'max_hvg': 8000, 'min_disp_mean': 0.3},
        'Drop-seq':     {'min_hvg': 500,  'max_hvg': 4000, 'min_disp_mean': 0.4},
        'CITE-seq':     {'min_hvg': 1000, 'max_hvg': 6000, 'min_disp_mean': 0.35},
        'Multiome':     {'min_hvg': 1000, 'max_hvg': 6000, 'min_disp_mean': 0.3},
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])

    metrics = {}
    warnings = []
    feedback = []

    if 'highly_variable' not in adata.var.columns:
        return {
            'metrics': {},
            'warnings': ['CRITICAL: highly_variable not in adata.var — HVG step may not have run'],
            'feedback': [],
            'success': False,
            'context_used': protocol
        }

    hvg_mask = adata.var['highly_variable']
    metrics['n_hvg'] = int(hvg_mask.sum())
    metrics['n_genes_total'] = int(adata.n_vars)
    metrics['pct_hvg'] = float(metrics['n_hvg'] / metrics['n_genes_total'] * 100)

    if 'dispersions_norm' in adata.var.columns:
        metrics['mean_dispersion_hvg'] = float(
            adata.var.loc[hvg_mask, 'dispersions_norm'].mean()
        )
    else:
        metrics['mean_dispersion_hvg'] = float('nan')

    # Check HVG count bounds
    if metrics['n_hvg'] < t['min_hvg']:
        warnings.append(
            f"TOO FEW HVGs: {metrics['n_hvg']} selected (min {t['min_hvg']} for {protocol}). "
            "Consider increasing n_top_genes or decreasing min_disp."
        )
        feedback.append('too_few_hvg')

    if metrics['n_hvg'] > t['max_hvg']:
        warnings.append(
            f"TOO MANY HVGs: {metrics['n_hvg']} selected (max {t['max_hvg']} for {protocol}). "
            "Consider decreasing n_top_genes."
        )
        feedback.append('too_many_hvg')

    # Check dispersion quality
    if (not float('nan') == metrics['mean_dispersion_hvg'] and
            metrics['mean_dispersion_hvg'] < t['min_disp_mean']):
        warnings.append(
            f"LOW DISPERSION: mean normalized dispersion of HVGs = "
            f"{metrics['mean_dispersion_hvg']:.3f} (threshold {t['min_disp_mean']}). "
            "Selected genes may be weakly variable."
        )
        feedback.append('low_variance_dominated')

    return {
        'metrics': metrics,
        'warnings': warnings,
        'feedback': feedback,
        'success': not bool(feedback),
        'context_used': protocol
    }
