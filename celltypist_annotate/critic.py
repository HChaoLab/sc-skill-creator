"""
celltypist_annotate critic — Evaluates cell type annotation quality

FEEDBACK KEYS (must match parameter_science_guide conditions in skill.yaml):
    high_unassigned       → p_thres too high or model mismatch
    low_confidence_labels → mean probability too low
    fragmented_voting     → min_prop too high; many mixed clusters unassigned
    over_homogenized      → too few cell types (possible over-annotation)
"""
import numpy as np


def critic_post_process(adata, context=None):
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')
    tissue   = context.get('tissue', 'immune')

    thresholds = {
        '10x Genomics': {'max_pct_unassigned': 20.0, 'min_mean_prob': 0.55, 'min_cell_types': 2},
        'Smart-seq2':   {'max_pct_unassigned': 15.0, 'min_mean_prob': 0.60, 'min_cell_types': 2},
        'Drop-seq':     {'max_pct_unassigned': 25.0, 'min_mean_prob': 0.50, 'min_cell_types': 2},
        'CITE-seq':     {'max_pct_unassigned': 20.0, 'min_mean_prob': 0.55, 'min_cell_types': 2},
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])

    metrics = {}
    warnings = []
    feedback = []

    # Find label column
    label_col = None
    for col in ['majority_voting', 'predicted_labels']:
        if col in adata.obs.columns:
            label_col = col
            break
    if label_col is None:
        return {
            'metrics': {},
            'warnings': ['CRITICAL: No CellTypist annotation columns found in adata.obs'],
            'feedback': [],
            'success': False,
            'context_used': protocol
        }

    labels = adata.obs[label_col]
    metrics['n_cells'] = int(len(labels))
    metrics['n_cell_types'] = int(labels.nunique())
    metrics['pct_unassigned'] = float((labels == 'Unassigned').sum() / len(labels) * 100)

    # Confidence score
    for conf_col in ['conf_score', 'max_prob']:
        if conf_col in adata.obs.columns:
            metrics['mean_max_probability'] = float(adata.obs[conf_col].mean())
            break
    else:
        metrics['mean_max_probability'] = float('nan')

    # Top cell types
    metrics['top3_celltypes'] = (
        labels[labels != 'Unassigned'].value_counts().head(3).index.tolist()
    )

    # Evaluate
    if metrics['pct_unassigned'] > t['max_pct_unassigned']:
        warnings.append(
            f"HIGH UNASSIGNED: {metrics['pct_unassigned']:.1f}% cells unassigned "
            f"(max {t['max_pct_unassigned']:.0f}% for {protocol}). "
            "Consider lowering p_thres or using a more appropriate model."
        )
        feedback.append('high_unassigned')

    if (not np.isnan(metrics['mean_max_probability']) and
            metrics['mean_max_probability'] < t['min_mean_prob']):
        warnings.append(
            f"LOW CONFIDENCE: mean probability {metrics['mean_max_probability']:.3f} "
            f"< threshold {t['min_mean_prob']}. "
            "Model may not match the tissue type. Check celltypist.models.models_description()."
        )
        feedback.append('low_confidence_labels')

    if metrics['n_cell_types'] < t['min_cell_types']:
        warnings.append(
            f"TOO FEW CELL TYPES: only {metrics['n_cell_types']} type(s) detected. "
            "Check if majority_voting is over-homogenizing results."
        )
        feedback.append('over_homogenized')

    return {
        'metrics': metrics,
        'warnings': warnings,
        'feedback': feedback,
        'success': not bool(feedback),
        'context_used': protocol
    }
