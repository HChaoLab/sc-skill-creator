"""
scanpy_rank_genes_groups critic — Evaluates marker gene detection quality

FEEDBACK KEYS (must match parameter_science_guide conditions in skill.yaml):
    too_few_markers      → n_genes too low; not enough markers for annotation
    low_specificity      → top scores are weak; markers not cluster-specific
    high_pct_background  → top markers expressed in many clusters; poor specificity
"""
import numpy as np


def critic_post_process(adata, context=None):
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')
    key = context.get('key_added', 'rank_genes_groups')

    thresholds = {
        '10x Genomics': {'min_top_score': 1.5,  'min_markers': 5,  'max_pct_background': 0.7},
        'Smart-seq2':   {'min_top_score': 2.0,  'min_markers': 10, 'max_pct_background': 0.6},
        'Drop-seq':     {'min_top_score': 1.0,  'min_markers': 5,  'max_pct_background': 0.7},
        'CITE-seq':     {'min_top_score': 1.5,  'min_markers': 5,  'max_pct_background': 0.7},
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])

    metrics = {}
    warnings = []
    feedback = []

    if key not in adata.uns:
        return {
            'metrics': {},
            'warnings': [f'CRITICAL: {key} not in adata.uns — run rank_genes_groups first'],
            'feedback': [],
            'success': False,
            'context_used': protocol
        }

    result = adata.uns[key]
    groups = list(result['names'].dtype.names)
    n_groups = len(groups)
    metrics['n_groups'] = n_groups

    top1_scores = []
    genes_per_group = []
    for g in groups:
        scores = result['scores'][g]
        if len(scores) > 0:
            top1_scores.append(float(scores[0]))
        genes_per_group.append(len(result['names'][g]))

    metrics['mean_top1_score'] = float(np.mean(top1_scores)) if top1_scores else 0.0
    metrics['min_top1_score'] = float(np.min(top1_scores)) if top1_scores else 0.0
    metrics['mean_genes_per_group'] = float(np.mean(genes_per_group))

    # Check marker quality
    if metrics['mean_genes_per_group'] < t['min_markers']:
        warnings.append(
            f"TOO FEW MARKERS: mean {metrics['mean_genes_per_group']:.0f} genes per group "
            f"(min {t['min_markers']}). Increase n_genes parameter."
        )
        feedback.append('too_few_markers')

    if metrics['mean_top1_score'] < t['min_top_score']:
        warnings.append(
            f"LOW MARKER SCORES: mean top-1 score = {metrics['mean_top1_score']:.2f} "
            f"(threshold {t['min_top_score']}). Markers may not be cluster-specific. "
            "Check clustering quality or switch method to 'wilcoxon'."
        )
        feedback.append('low_specificity')

    # Check for background expression if pts available
    if 'pts_rest' in result:
        high_bg_count = 0
        for g in groups:
            pts_rest = result['pts_rest'][g]
            if len(pts_rest) > 0 and float(pts_rest[0]) > t['max_pct_background']:
                high_bg_count += 1
        if high_bg_count > n_groups // 2:
            warnings.append(
                f"HIGH BACKGROUND: top markers in {high_bg_count}/{n_groups} groups "
                f"are expressed in >{t['max_pct_background']:.0%} of other cells. "
                "Markers lack cluster specificity."
            )
            feedback.append('high_pct_background')

    return {
        'metrics': metrics,
        'warnings': warnings,
        'feedback': feedback,
        'success': not bool(feedback),
        'context_used': protocol
    }
