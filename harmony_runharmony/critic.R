#' harmony_runharmony critic — Evaluates Harmony batch correction quality
#'
#' FEEDBACK KEYS (match parameter_science_guide in skill.yaml):
#'   under_corrected         → batches still separate; theta too low
#'   over_corrected          → biological structure lost; theta too high
#'   convergence_not_reached → max.iter.harmony too low

critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol  <- if (!is.null(context$protocol))  context$protocol  else "10x Genomics"
    n_batches <- if (!is.null(context$n_batches))  context$n_batches else NA_integer_

    thresholds <- list(
        "10x Genomics" = list(min_mixing = 0.3, max_mixing = 0.95),
        "Smart-seq2"   = list(min_mixing = 0.2, max_mixing = 0.90),
        "Drop-seq"     = list(min_mixing = 0.3, max_mixing = 0.95),
        "CITE-seq"     = list(min_mixing = 0.25, max_mixing = 0.90)
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]

    metrics  <- list()
    warnings <- character(0)
    feedback <- character(0)

    # Check harmony reduction exists
    if (!"harmony" %in% names(obj@reductions)) {
        return(list(
            metrics      = list(),
            warnings     = "CRITICAL: harmony reduction not found — run RunHarmony() first",
            feedback     = character(0),
            success      = FALSE,
            context_used = protocol
        ))
    }

    # Extract metrics from analysis_history
    batch_mixing_score <- NA_real_
    history <- obj@misc$analysis_history
    if (!is.null(history)) {
        for (entry in rev(history)) {
            if (!is.null(entry$skill_id) && entry$skill_id == "harmony_runharmony") {
                batch_mixing_score  <- entry$metrics$batch_mixing_score
                metrics$n_batches   <- entry$metrics$n_batches
                metrics$n_cells     <- entry$metrics$n_cells
                metrics$n_harmony_dims <- entry$metrics$n_harmony_dims
                break
            }
        }
    }
    metrics$batch_mixing_score <- batch_mixing_score

    # Evaluate batch mixing
    if (!is.na(batch_mixing_score)) {
        if (batch_mixing_score < t$min_mixing) {
            warnings <- c(warnings, sprintf(
                "UNDER-CORRECTED: batch mixing score %.2f < threshold %.2f for %s. ",
                batch_mixing_score, t$min_mixing, protocol
            ))
            feedback <- c(feedback, "under_corrected")
        }
        if (batch_mixing_score > t$max_mixing) {
            warnings <- c(warnings, sprintf(
                "POSSIBLE OVER-CORRECTION: batch mixing score %.2f > %.2f. ",
                batch_mixing_score, t$max_mixing
            ))
            feedback <- c(feedback, "over_corrected")
        }
    }

    # Check if Harmony UMAP still shows batch separation
    # (proxy: check if leiden/louvain clusters are batch-dominated)
    cluster_cols <- intersect(c("leiden", "louvain", "seurat_clusters"), colnames(obj@meta.data))
    if (length(cluster_cols) > 0 && !is.null(history)) {
        # Find group.by.vars from history
        for (entry in rev(history)) {
            if (!is.null(entry$skill_id) && entry$skill_id == "harmony_runharmony") {
                batch_col <- entry$params$group.by.vars[1]
                if (!is.null(batch_col) && batch_col %in% colnames(obj@meta.data)) {
                    cluster_tab <- table(
                        obj@meta.data[[cluster_cols[1]]],
                        obj@meta.data[[batch_col]]
                    )
                    # Compute entropy-based mixing per cluster
                    if (ncol(cluster_tab) > 1) {
                        row_entropies <- apply(cluster_tab, 1, function(row) {
                            p <- row / sum(row)
                            p <- p[p > 0]
                            -sum(p * log(p))
                        })
                        max_entropy <- log(ncol(cluster_tab))
                        mean_rel_entropy <- mean(row_entropies) / max_entropy
                        metrics$mean_cluster_batch_entropy <- mean_rel_entropy

                        if (mean_rel_entropy < 0.3) {
                            warnings <- c(warnings, sprintf(
                                "BATCH-DOMINATED CLUSTERS: mean relative entropy %.2f. ",
                                mean_rel_entropy
                            ))
                            if (!"under_corrected" %in% feedback) {
                                feedback <- c(feedback, "under_corrected")
                            }
                        }
                    }
                }
                break
            }
        }
    }

    list(
        metrics      = metrics,
        warnings     = warnings,
        feedback     = feedback,
        success      = length(feedback) == 0L,
        context_used = protocol
    )
}
