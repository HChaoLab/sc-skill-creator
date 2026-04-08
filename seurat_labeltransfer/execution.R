#' seurat_labeltransfer — Seurat Label Transfer
#'
#' PURPOSE:
#'   Transfers cell type annotations from an annotated reference Seurat object
#'   to the query dataset using CCA anchor-based integration.
#'
#' INPUTS:
#'   input_data  : character (RDS path) | Seurat object (query)
#'   params_list : list (optional, overrides defaults)
#'                 REQUIRED: params_list$reference_path = path to reference RDS
#'
#' DEFAULTS:
#'   ref_label_col = "cell_type", dims = 1:30,
#'   normalization_method = "LogNormalize",
#'   k_filter = 200, k_score = 30, k_weight = 100
#'
#' OUTPUTS:
#'   Query Seurat with predicted labels in @meta.data and @misc$analysis_history

library(Seurat)

run_labeltransfer <- function(input_data, params_list = NULL) {
    # Memory-First
    if (is.character(input_data)) {
        query <- readRDS(input_data)
    } else {
        query <- input_data
    }

    # Parameter Safety
    default_params <- list(
        reference_path       = NULL,
        ref_label_col        = "cell_type",
        dims                 = 1:30,
        normalization_method = "LogNormalize",
        k_filter             = 200L,
        k_score              = 30L,
        k_weight             = 100L,
        prediction_col       = "predicted_celltype"
    )
    current_params <- modifyList(
        default_params,
        if (!is.null(params_list)) params_list else list()
    )

    # Validate reference
    if (is.null(current_params$reference_path)) {
        stop("reference_path is required. Provide path to reference Seurat RDS file in params_list.")
    }
    reference <- readRDS(current_params$reference_path)

    # Validate label column
    if (!current_params$ref_label_col %in% colnames(reference@meta.data)) {
        available <- colnames(reference@meta.data)
        stop(paste0(
            "Label column '", current_params$ref_label_col, "' not found in reference@meta.data. ",
            "Available columns: ", paste(available, collapse = ", ")
        ))
    }

    # Find transfer anchors
    anchors <- Seurat::FindTransferAnchors(
        reference            = reference,
        query                = query,
        dims                 = current_params$dims,
        normalization.method = current_params$normalization_method,
        k.filter             = as.integer(current_params$k_filter),
        k.score              = as.integer(current_params$k_score)
    )

    # Transfer labels
    predictions <- Seurat::TransferData(
        anchorset = anchors,
        refdata   = reference@meta.data[[current_params$ref_label_col]],
        dims      = current_params$dims,
        k.weight  = as.integer(current_params$k_weight)
    )

    # Add to query
    pred_col    <- current_params$prediction_col
    score_col   <- paste0(pred_col, ".score")
    query[[pred_col]]  <- predictions$predicted.id
    query[[score_col]] <- predictions$prediction.score.max

    # Collect metrics
    labels           <- query@meta.data[[pred_col]]
    scores           <- query@meta.data[[score_col]]
    n_predicted_types <- length(unique(labels))
    mean_pred_score  <- mean(scores, na.rm = TRUE)
    pct_high_conf    <- mean(scores >= 0.75, na.rm = TRUE) * 100
    top3             <- names(sort(table(labels), decreasing = TRUE))[1:min(3L, n_predicted_types)]

    # Metadata Footprinting
    if (is.null(query@misc$analysis_history)) query@misc$analysis_history <- list()
    query@misc$analysis_history <- c(query@misc$analysis_history, list(list(
        skill_id  = "seurat_labeltransfer",
        params    = list(
            ref_label_col        = current_params$ref_label_col,
            dims                 = as.integer(current_params$dims),
            normalization_method = current_params$normalization_method,
            k_filter             = as.integer(current_params$k_filter),
            k_score              = as.integer(current_params$k_score),
            k_weight             = as.integer(current_params$k_weight)
        ),
        metrics   = list(
            n_predicted_types  = as.integer(n_predicted_types),
            mean_prediction_score = mean_pred_score,
            pct_high_confidence   = pct_high_conf,
            top3_predicted        = top3,
            n_cells               = as.integer(ncol(query)),
            n_anchors             = nrow(anchors@anchors)
        ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )))

    if (is.character(input_data)) saveRDS(query, input_data)
    return(query)
}
