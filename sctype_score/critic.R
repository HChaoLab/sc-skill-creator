#' sctype_score critic — Evaluates ScType cell type assignment quality
#'
#' FEEDBACK KEYS (match parameter_science_guide in skill.yaml):
#'   high_unknown  → many clusters labelled Unknown; tissue_type mismatch
#'   low_scores    → mean scores low; wrong assay or tissue
#'   too_few_types → fewer cell types than clusters; over-collapsed annotation

critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol    <- if (!is.null(context$protocol))    context$protocol    else "10x Genomics"
    tissue_type <- if (!is.null(context$tissue_type)) context$tissue_type else "Immune system"

    thresholds <- list(
        "10x Genomics" = list(max_pct_unknown = 30.0, min_score = 0.0, min_types = 2L),
        "Smart-seq2"   = list(max_pct_unknown = 20.0, min_score = 0.0, min_types = 2L),
        "Drop-seq"     = list(max_pct_unknown = 35.0, min_score = 0.0, min_types = 2L),
        "CITE-seq"     = list(max_pct_unknown = 30.0, min_score = 0.0, min_types = 2L)
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]

    metrics  <- list()
    warnings <- character(0)
    feedback <- character(0)

    # Find ScType column
    sctype_cols <- grep("sctype_classification", colnames(obj@meta.data), value = TRUE)
    sctype_cols <- sctype_cols[!grepl("score", sctype_cols)]

    if (length(sctype_cols) == 0L) {
        return(list(
            metrics      = list(),
            warnings     = "CRITICAL: No ScType classification columns found — run sctype_score first",
            feedback     = character(0),
            success      = FALSE,
            context_used = protocol
        ))
    }

    label_col  <- sctype_cols[1]
    score_col  <- paste0(label_col, ".score")
    labels     <- obj@meta.data[[label_col]]
    has_scores <- score_col %in% colnames(obj@meta.data)

    metrics$n_cells          <- as.integer(ncol(obj))
    metrics$n_clusters       <- as.integer(length(unique(
        obj@meta.data[["seurat_clusters"]])))
    metrics$n_assigned_types <- as.integer(length(unique(labels[labels != "Unknown"])))
    metrics$pct_unknown      <- mean(labels == "Unknown") * 100

    if (has_scores) {
        valid_scores          <- obj@meta.data[[score_col]][labels != "Unknown"]
        metrics$mean_top_score <- mean(valid_scores, na.rm = TRUE)
    } else {
        metrics$mean_top_score <- NA_real_
    }

    metrics$top3_types <- names(
        sort(table(labels[labels != "Unknown"]), decreasing = TRUE)
    )[1:min(3L, metrics$n_assigned_types)]

    if (metrics$pct_unknown > t$max_pct_unknown) {
        warnings <- c(warnings, sprintf(
            "HIGH UNKNOWN RATE: %.1f%% cells classified as Unknown (max %.0f%% for %s). "
            "Check tissue_type or provide a custom marker database.",
            metrics$pct_unknown, t$max_pct_unknown, protocol
        ))
        feedback <- c(feedback, "high_unknown")
    }

    if (!is.na(metrics$mean_top_score) && metrics$mean_top_score < t$min_score) {
        warnings <- c(warnings, sprintf(
            "LOW SCTYPE SCORES: mean score = %.4f. "
            "Ensure scaled_assay is correct (scale.data vs SCT).",
            metrics$mean_top_score
        ))
        feedback <- c(feedback, "low_scores")
    }

    if (metrics$n_assigned_types < t$min_types) {
        warnings <- c(warnings, sprintf(
            "TOO FEW CELL TYPES: %d type(s) assigned. "
            "Try a more specific tissue_type or custom marker database.",
            metrics$n_assigned_types
        ))
        feedback <- c(feedback, "too_few_types")
    }

    list(
        metrics      = metrics,
        warnings     = warnings,
        feedback     = feedback,
        success      = length(feedback) == 0L,
        context_used = protocol
    )
}
