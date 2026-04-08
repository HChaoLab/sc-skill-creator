"""
scanpy_umap critic — Evaluates UMAP embedding quality

FEEDBACK KEYS (must match parameter_science_guide conditions in skill.yaml):
    overlapping_clusters → min_dist too high; clusters not well-separated
    too_compressed       → embedding range very small; try increasing spread
    clusters_too_close   → spread too small; clusters visually merged
    clusters_too_spread  → spread too large; embedding overly dispersed
"""
import numpy as np


def critic_post_process(adata, context=None):
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')

    thresholds = {
        '10x Genomics': {'min_range': 5.0,  'max_range': 60.0},
        'Smart-seq2':   {'min_range': 3.0,  'max_range': 50.0},
        'Drop-seq':     {'min_range': 5.0,  'max_range': 55.0},
        'CITE-seq':     {'min_range': 4.0,  'max_range': 60.0},
        'Multiome':     {'min_range': 4.0,  'max_range': 60.0},
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])

    metrics = {}
    warnings = []
    feedback = []

    if 'X_umap' not in adata.obsm:
        return {
            'metrics': {},
            'warnings': ['CRITICAL: X_umap not found — run sc.tl.umap() first'],
            'feedback': [],
            'success': False,
            'context_used': protocol
        }

    umap = adata.obsm['X_umap']
    range_x = float(umap[:, 0].max() - umap[:, 0].min())
    range_y = float(umap[:, 1].max() - umap[:, 1].min())
    metrics['umap_range_x'] = range_x
    metrics['umap_range_y'] = range_y
    metrics['umap_center_x'] = float(umap[:, 0].mean())
    metrics['umap_center_y'] = float(umap[:, 1].mean())

    if range_x < t['min_range'] or range_y < t['min_range']:
        warnings.append(
            f"COMPRESSED EMBEDDING: UMAP range ({range_x:.1f} x {range_y:.1f}) is very small "
            f"(min {t['min_range']} for {protocol}). Try increasing spread or decreasing min_dist."
        )
        feedback.append('too_compressed')

    if range_x > t['max_range'] or range_y > t['max_range']:
        warnings.append(
            f"OVER-DISPERSED EMBEDDING: UMAP range ({range_x:.1f} x {range_y:.1f}) is very large "
            f"(max {t['max_range']}). Consider decreasing spread."
        )
        feedback.append('clusters_too_spread')

    # Check cluster separation if leiden/louvain present
    for cluster_key in ['leiden', 'louvain']:
        if cluster_key in adata.obs.columns:
            labels = adata.obs[cluster_key]
            if labels.nunique() > 1:
                # Compute mean pairwise centroid distance
                centroids = np.array([
                    umap[labels == c].mean(axis=0) for c in labels.unique()
                ])
                if len(centroids) > 1:
                    from itertools import combinations
                    dists = [
                        np.linalg.norm(centroids[i] - centroids[j])
                        for i, j in combinations(range(len(centroids)), 2)
                    ]
                    metrics['mean_centroid_distance'] = float(np.mean(dists))
                    if metrics['mean_centroid_distance'] < 1.5:
                        warnings.append(
                            f"OVERLAPPING CLUSTERS: mean centroid distance "
                            f"{metrics['mean_centroid_distance']:.2f} is very small. "
                            "Try decreasing min_dist."
                        )
                        feedback.append('overlapping_clusters')
            break

    return {
        'metrics': metrics,
        'warnings': warnings,
        'feedback': feedback,
        'success': not bool(feedback),
        'context_used': protocol
    }
