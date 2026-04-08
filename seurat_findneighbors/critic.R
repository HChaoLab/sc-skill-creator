#' seurat_findneighbors critic — Evaluates SNN graph quality
#'
#' FEEDBACK KEYS (match parameter_science_guide in skill.yaml):
#'   missing_rare_populations → dims too few; rare variation excluded
#'   noisy_clusters           → dims too many; noise in distance calc
#'   fragmented_clusters      → k.param too low; graph under-connected
#'   over_merged_populations  → k.param too high; populations merged

critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol <- if (!is.null(context$protocol)) context$protocol else "10x Genomics"

    thresholds <- list(
        "10x Genomics" = list(min_dims = 5L,  max_dims = 50L, min_k = 5L,  max_k = 50L),
        "Smart-seq2"   = list(min_dims = 10L, max_dims = 50L, min_k = 10L, max_k = 50L),
        "Drop-seq"     = list(min_dims = 5L,  max_dims = 40L, min_k = 5L,  max_k = 40L),
        "CITE-seq"     = list(min_dims = 10L, max_dims = 50L, min_k = 10L, max_k = 50L)
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]

    metrics  <- list()
    warnings <- character(0)
    feedback <- character(0)

    # Check graph exists
    graph_names <- names(obj@graphs)
    snn_found <- any(grepl("snn", tolower(graph_names)))
    if (!snn_found) {
        return(list(
            metrics      = list(),
            warnings     = "CRITICAL: No SNN graph found — run FindNeighbors() first",
            feedback     = character(0),
            success      = FALSE,
            context_used = protocol
        ))
    }

    # Extract params from analysis_history if available
    history <- obj@misc$analysis_history
    n_dims_used  <- NA_integer_
    k_param_used <- NA_integer_

    if (!is.null(history) && length(history) > 0) {
        for (entry in rev(history)) {
            if (!is.null(entry$skill_id) && entry$skill_id == "seurat_findneighbors") {
                n_dims_used  <- entry$metrics$n_dims_used
                k_param_used <- entry$metrics$k_param_used
                break
            }
        }
    }

    metrics$n_cells      <- as.integer(ncol(obj))
    metrics$n_dims_used  <- n_dims_used
    metrics$k_param_used <- k_param_used

    if (!is.na(n_dims_used)) {
        if (n_dims_used < t$min_dims) {
            warnings <- c(warnings, sprintf(
                "TOO FEW DIMS: %d dims used (min %d for %s). Rare populations may be missed.",
                n_dims_used, t$min_dims, protocol
            ))
            feedback <- c(feedback, "missing_rare_populations")
        }
        if (n_dims_used > t$max_dims) {
            warnings <- c(warnings, sprintf(
                "TOO MANY DIMS: %d dims used (max %d for %s). Noise may dominate graph.",
                n_dims_used, t$max_dims, protocol
            ))
            feedback <- c(feedback, "noisy_clusters")
        }
    }

    if (!is.na(k_param_used)) {
        if (k_param_used < t$min_k) {
            warnings <- c(warnings, sprintf(
                "LOW k.param: %d (min %d for %s). Graph may be fragmented.",
                k_param_used, t$min_k, protocol
            ))
            feedback <- c(feedback, "fragmented_clusters")
        }
        if (k_param_used > t$max_k) {
            warnings <- c(warnings, sprintf(
                "HIGH k.param: %d (max %d for %s). Distinct populations may be merged.",
                k_param_used, t$max_k, protocol
            ))
            feedback <- c(feedback, "over_merged_populations")
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
