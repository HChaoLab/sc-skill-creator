#' harmony_runharmony — Harmony Batch Correction
#'
#' PURPOSE:
#'   Corrects batch effects by aligning PCA embeddings across batches using
#'   iterative clustering. Only modifies the 'harmony' reduction — expression
#'   values are untouched. Use 'harmony' reduction in FindNeighbors/RunUMAP.
#'
#' INPUTS:
#'   input_data  : character (RDS path) | Seurat object
#'   params_list : list (optional, overrides defaults)
#'
#' DEFAULTS:
#'   group.by.vars = "orig.ident", theta = 2.0, lambda = 1.0,
#'   sigma = 0.1, nclust = NULL, max.iter.harmony = 10, dims.use = NULL
#'
#' OUTPUTS:
#'   Seurat object with @reductions$harmony and @misc$analysis_history updated
#'   NOTE: Use reduction = "harmony" in subsequent FindNeighbors() and RunUMAP()

library(harmony)

run_runharmony <- function(input_data, params_list = NULL) {
    # Memory-First
    if (is.character(input_data)) {
        obj <- readRDS(input_data)
    } else {
        obj <- input_data
    }

    # Parameter Safety
    default_params <- list(
        group.by.vars     = "orig.ident",
        theta             = 2.0,
        lambda            = 1.0,
        sigma             = 0.1,
        nclust            = NULL,
        max.iter.harmony  = 10L,
        dims.use          = NULL
    )
    current_params <- modifyList(
        default_params,
        if (!is.null(params_list)) params_list else list()
    )

    # Validate group.by.vars
    for (col in current_params$group.by.vars) {
        if (!col %in% colnames(obj@meta.data)) {
            stop(paste0(
                "Column '", col, "' not found in obj@meta.data. ",
                "Available: ", paste(colnames(obj@meta.data), collapse = ", ")
            ))
        }
    }

    # Check PCA exists
    if (!"pca" %in% names(obj@reductions)) {
        stop("PCA reduction not found. Run RunPCA() (or SCTransform + RunPCA) before RunHarmony().")
    }

    # Build args
    harmony_args <- list(
        object           = obj,
        group.by.vars    = current_params$group.by.vars,
        theta            = current_params$theta,
        lambda           = current_params$lambda,
        sigma            = current_params$sigma,
        max.iter.harmony = as.integer(current_params$max.iter.harmony),
        verbose          = FALSE
    )
    if (!is.null(current_params$nclust)) {
        harmony_args$nclust <- as.integer(current_params$nclust)
    }
    if (!is.null(current_params$dims.use)) {
        harmony_args$dims.use <- as.integer(current_params$dims.use)
    }

    obj <- do.call(harmony::RunHarmony, harmony_args)

    # Collect metrics
    harmony_embed <- obj@reductions[["harmony"]]@cell.embeddings
    n_harmony_dims <- ncol(harmony_embed)
    n_cells <- nrow(harmony_embed)

    # Batch mixing score: mean silhouette-like metric across batches
    batch_labels <- obj@meta.data[[current_params$group.by.vars[1]]]
    n_batches    <- length(unique(batch_labels))

    # Simple mixing proxy: SD of batch centroids (lower = better mixed)
    batch_mixing_score <- NA_real_
    if (n_batches > 1 && n_batches <= 50) {
        centroids <- sapply(unique(batch_labels), function(b) {
            colMeans(harmony_embed[batch_labels == b, 1:min(2, n_harmony_dims), drop = FALSE])
        })
        if (is.matrix(centroids)) {
            centroid_sd <- mean(apply(centroids, 1, sd))
            # Normalize by overall embedding SD
            overall_sd <- mean(apply(harmony_embed[, 1:min(2, n_harmony_dims), drop = FALSE], 2, sd))
            batch_mixing_score <- if (overall_sd > 0) 1 - centroid_sd / overall_sd else NA_real_
            batch_mixing_score <- max(0, min(1, batch_mixing_score))
        }
    }

    # Metadata Footprinting
    if (is.null(obj@misc$analysis_history)) obj@misc$analysis_history <- list()
    obj@misc$analysis_history <- c(obj@misc$analysis_history, list(list(
        skill_id  = "harmony_runharmony",
        params    = list(
            group.by.vars    = current_params$group.by.vars,
            theta            = current_params$theta,
            lambda           = current_params$lambda,
            sigma            = current_params$sigma,
            max.iter.harmony = as.integer(current_params$max.iter.harmony)
        ),
        metrics   = list(
            n_harmony_clusters = NA_integer_,   # Harmony doesn't expose cluster count directly
            n_harmony_dims     = as.integer(n_harmony_dims),
            n_batches          = as.integer(n_batches),
            n_cells            = as.integer(n_cells),
            batch_mixing_score = batch_mixing_score
        ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )))

    if (is.character(input_data)) saveRDS(obj, input_data)
    return(obj)
}
