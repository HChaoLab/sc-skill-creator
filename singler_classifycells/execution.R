#' singler_classifycells — SingleR Reference-Based Cell Type Annotation
#'
#' PURPOSE:
#'   Annotates cell types by correlating expression profiles against curated
#'   reference datasets (celldex). Gold-standard automated annotation method.
#'   Accepts Seurat or SingleCellExperiment objects.
#'
#' INPUTS:
#'   input_data  : character (RDS path) | Seurat | SingleCellExperiment
#'   params_list : list (optional, overrides defaults)
#'
#' DEFAULTS:
#'   ref = "HumanPrimaryCellAtlasData", labels_field = "label.main",
#'   assay_type = "logcounts", de_method = "classic",
#'   aggr_ref = FALSE, prune = TRUE
#'
#' OUTPUTS:
#'   Object with SingleR labels in colData/meta.data and @misc$analysis_history

library(SingleR)
library(celldex)

run_classifycells <- function(input_data, params_list = NULL) {
    # Memory-First
    if (is.character(input_data)) {
        obj <- readRDS(input_data)
    } else {
        obj <- input_data
    }

    # Parameter Safety
    default_params <- list(
        ref          = "HumanPrimaryCellAtlasData",
        labels_field = "label.main",
        assay_type   = "logcounts",
        de_method    = "classic",
        aggr_ref     = FALSE,
        prune        = TRUE
    )
    current_params <- modifyList(
        default_params,
        if (!is.null(params_list)) params_list else list()
    )

    # Convert Seurat → SCE if needed
    is_seurat <- inherits(obj, "Seurat")
    if (is_seurat) {
        if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat required for conversion")
        sce <- Seurat::as.SingleCellExperiment(obj)
        # Seurat's SCE uses 'logcounts' from the data slot
    } else {
        sce <- obj
    }

    # Load reference
    ref_fn_name <- current_params$ref
    ref_fn <- tryCatch(
        get(ref_fn_name, envir = asNamespace("celldex")),
        error = function(e) stop(paste0(
            "Reference '", ref_fn_name, "' not found in celldex. ",
            "Available: HumanPrimaryCellAtlasData, BlueprintEncodeData, ",
            "ImmGenData, MouseRNAseqData, DatabaseImmuneCellExpressionData, etc."
        ))
    )
    ref_data <- ref_fn()

    # Run SingleR
    singler_result <- SingleR::SingleR(
        test        = sce,
        ref         = ref_data,
        labels      = ref_data[[current_params$labels_field]],
        de.method   = current_params$de_method,
        aggr.ref    = current_params$aggr_ref,
        assay.type.test = current_params$assay_type
    )

    # Pruning
    if (current_params$prune) {
        singler_result <- SingleR::pruneScores(singler_result)
    }

    # Add labels back to object
    final_labels <- if (current_params$prune) {
        singler_result$pruned.labels
    } else {
        singler_result$labels
    }

    label_col <- paste0("SingleR.", gsub("\\.", "_", current_params$labels_field))

    if (is_seurat) {
        obj@meta.data[[label_col]] <- final_labels
        obj@meta.data[["SingleR.score"]] <- singler_result$scores[
            cbind(seq_len(nrow(singler_result$scores)),
                  match(singler_result$labels, colnames(singler_result$scores)))
        ]
        obj@meta.data[["SingleR.delta.next"]] <- singler_result$delta.next
    } else {
        SummarizedExperiment::colData(obj)[[label_col]] <- final_labels
        SummarizedExperiment::colData(obj)[["SingleR.delta.next"]] <- singler_result$delta.next
    }

    # Collect metrics
    n_cell_types <- length(unique(na.omit(final_labels)))
    pct_pruned   <- if (current_params$prune) {
        mean(is.na(singler_result$pruned.labels)) * 100
    } else {
        0.0
    }
    mean_delta_next <- mean(singler_result$delta.next, na.rm = TRUE)
    top3 <- names(sort(table(na.omit(final_labels)), decreasing = TRUE))[1:min(3, n_cell_types)]

    # Metadata Footprinting — store in misc for Seurat, metadata() for SCE
    footprint <- list(
        skill_id  = "singler_classifycells",
        params    = list(
            ref          = current_params$ref,
            labels_field = current_params$labels_field,
            de_method    = current_params$de_method,
            aggr_ref     = current_params$aggr_ref,
            prune        = current_params$prune
        ),
        metrics   = list(
            n_cell_types    = as.integer(n_cell_types),
            pct_pruned      = pct_pruned,
            mean_delta_next = mean_delta_next,
            top3_celltypes  = top3,
            n_cells         = as.integer(ncol(sce))
        ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )

    if (is_seurat) {
        if (is.null(obj@misc$analysis_history)) obj@misc$analysis_history <- list()
        obj@misc$analysis_history <- c(obj@misc$analysis_history, list(footprint))
    } else {
        hist <- S4Vectors::metadata(obj)$analysis_history
        if (is.null(hist)) hist <- list()
        S4Vectors::metadata(obj)$analysis_history <- c(hist, list(footprint))
    }

    if (is.character(input_data)) saveRDS(obj, input_data)
    return(obj)
}
