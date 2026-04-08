#' seurat_labeltransfer critic — Evaluates label transfer quality
#'
#' FEEDBACK KEYS (match parameter_science_guide in skill.yaml):
#'   low_prediction_score  → mean score too low; try more dims or different reference
#'   noisy_transfer        → high score variance; too many dims
#'   fragmented_predictions → predictions inconsistent within clusters
#'   over_smoothed         → too few types detected; k_weight too high

critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol <- if (!is.null(context$protocol)) context$protocol else "10x Genomics"

    thresholds <- list(
        "10x Genomics" = list(min_pred_score = 0.5, min_high_conf = 40.0, min_types = 2L),
        "Smart-seq2"   = list(min_pred_score = 0.6, min_high_conf = 50.0, min_types = 2L),
        "Drop-seq"     = list(min_pred_score = 0.45, min_high_conf = 35.0, min_types = 2L),
        "CITE-seq"     = list(min_pred_score = 0.55, min_high_conf = 45.0, min_types = 2L)
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]

    metrics  <- list()
    warnings <- character(0)
    feedback <- character(0)

    # Find prediction column
    pred_cols <- grep("predicted_celltype|predicted\\.id", colnames(obj@meta.data), value = TRUE)
    score_cols <- grep("\\.score$", colnames(obj@meta.data), value = TRUE)

    if (length(pred_cols) == 0L) {
        return(list(
            metrics      = list(),
            warnings     = "CRITICAL: No label transfer prediction columns found — run seurat_labeltransfer first",
            feedback     = character(0),
            success      = FALSE,
            context_used = protocol
        ))
    }

    labels <- obj@meta.data[[pred_cols[1]]]
    metrics$n_cells           <- as.integer(ncol(obj))
    metrics$n_predicted_types <- as.integer(length(unique(labels)))
    metrics$top3_predicted    <- names(sort(table(labels), decreasing = TRUE))[1:min(3L, metrics$n_predicted_types)]

    if (length(score_cols) > 0) {
        scores <- obj@meta.data[[score_cols[1]]]
        metrics$mean_prediction_score <- mean(scores, na.rm = TRUE)
        metrics$pct_high_confidence   <- mean(scores >= 0.75, na.rm = TRUE) * 100
        metrics$sd_prediction_score   <- sd(scores, na.rm = TRUE)

        if (metrics$mean_prediction_score < t$min_pred_score) {
            warnings <- c(warnings, sprintf(
                "LOW PREDICTION SCORE: mean = %.3f (threshold %.2f for %s). "
                "Reference may not match tissue. Try increasing dims.",
                metrics$mean_prediction_score, t$min_pred_score, protocol
            ))
            feedback <- c(feedback, "low_prediction_score")
        }

        if (metrics$pct_high_confidence < t$min_high_conf) {
            warnings <- c(warnings, sprintf(
                "FEW HIGH-CONFIDENCE PREDICTIONS: %.1f%% cells with score ≥ 0.75 "
                "(threshold %.0f%%). Reference may be mismatched.",
                metrics$pct_high_confidence, t$min_high_conf
            ))
            if (!"low_prediction_score" %in% feedback) {
                feedback <- c(feedback, "low_prediction_score")
            }
        }

        # Check prediction consistency within clusters
        cluster_cols <- intersect(c("seurat_clusters", "leiden", "louvain"), colnames(obj@meta.data))
        if (length(cluster_cols) > 0) {
            cluster_labels <- table(obj@meta.data[[cluster_cols[1]]], labels)
            cluster_entropy <- apply(cluster_labels, 1, function(row) {
                p <- row / sum(row)
                p <- p[p > 0]
                -sum(p * log(p))
            })
            if (mean(cluster_entropy) > 1.5) {
                warnings <- c(warnings,
                    "FRAGMENTED PREDICTIONS: High label entropy within clusters. "
                    "Consider increasing k_weight for smoother predictions."
                )
                feedback <- c(feedback, "fragmented_predictions")
            }
        }
    }

    if (metrics$n_predicted_types < t$min_types) {
        warnings <- c(warnings, sprintf(
            "TOO FEW CELL TYPES: %d type(s) predicted. k_weight may be too high.",
            metrics$n_predicted_types
        ))
        feedback <- c(feedback, "over_smoothed")
    }

    list(
        metrics      = metrics,
        warnings     = warnings,
        feedback     = feedback,
        success      = length(feedback) == 0L,
        context_used = protocol
    )
}
