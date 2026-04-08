#' singler_classifycells critic — Evaluates SingleR annotation quality
#'
#' FEEDBACK KEYS (match parameter_science_guide in skill.yaml):
#'   too_few_cell_types  → labels_field too coarse; try label.fine
#'   too_many_pruned     → high ambiguity; switch to label.main
#'   high_na_rate        → prune=TRUE removing too many cells
#'   low_confidence_labels → mean delta.next too low

critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol <- if (!is.null(context$protocol)) context$protocol else "10x Genomics"

    thresholds <- list(
        "10x Genomics" = list(max_pct_pruned = 20.0, min_delta_next = 0.05, min_cell_types = 2L),
        "Smart-seq2"   = list(max_pct_pruned = 15.0, min_delta_next = 0.08, min_cell_types = 2L),
        "Drop-seq"     = list(max_pct_pruned = 25.0, min_delta_next = 0.04, min_cell_types = 2L),
        "CITE-seq"     = list(max_pct_pruned = 20.0, min_delta_next = 0.05, min_cell_types = 2L)
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]

    metrics  <- list()
    warnings <- character(0)
    feedback <- character(0)

    # Find SingleR label column
    is_seurat <- inherits(obj, "Seurat")
    meta <- if (is_seurat) obj@meta.data else as.data.frame(SummarizedExperiment::colData(obj))
    singler_cols <- grep("^SingleR\\.", colnames(meta), value = TRUE)

    if (length(singler_cols) == 0L) {
        return(list(
            metrics      = list(),
            warnings     = "CRITICAL: No SingleR annotation columns found — run SingleR first",
            feedback     = character(0),
            success      = FALSE,
            context_used = protocol
        ))
    }

    label_col <- singler_cols[!grepl("score|delta", singler_cols)][1]
    labels <- meta[[label_col]]

    metrics$n_cells      <- as.integer(nrow(meta))
    metrics$n_cell_types <- as.integer(length(unique(na.omit(labels))))
    metrics$pct_pruned   <- mean(is.na(labels)) * 100

    if ("SingleR.delta.next" %in% colnames(meta)) {
        metrics$mean_delta_next <- mean(meta$SingleR.delta.next, na.rm = TRUE)
    } else {
        metrics$mean_delta_next <- NA_real_
    }

    metrics$top3_celltypes <- names(
        sort(table(na.omit(labels)), decreasing = TRUE)
    )[1:min(3L, metrics$n_cell_types)]

    if (metrics$pct_pruned > t$max_pct_pruned) {
        warnings <- c(warnings, sprintf(
            "HIGH PRUNING RATE: %.1f%% cells pruned (max %.0f%% for %s). "
            "Consider switching to 'label.main' for broader categories.",
            metrics$pct_pruned, t$max_pct_pruned, protocol
        ))
        feedback <- c(feedback, "too_many_pruned")
    }

    if (!is.na(metrics$mean_delta_next) && metrics$mean_delta_next < t$min_delta_next) {
        warnings <- c(warnings, sprintf(
            "LOW CONFIDENCE: mean delta.next = %.3f (threshold %.3f). "
            "Many cells are ambiguous between closely related types.",
            metrics$mean_delta_next, t$min_delta_next
        ))
        feedback <- c(feedback, "low_confidence_labels")
    }

    if (metrics$n_cell_types < t$min_cell_types) {
        warnings <- c(warnings, sprintf(
            "TOO FEW CELL TYPES: %d type(s) detected. "
            "Try 'label.fine' for more granular annotation.",
            metrics$n_cell_types
        ))
        feedback <- c(feedback, "too_few_cell_types")
    }

    list(
        metrics      = metrics,
        warnings     = warnings,
        feedback     = feedback,
        success      = length(feedback) == 0L,
        context_used = protocol
    )
}
