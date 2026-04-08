#' seurat_findneighbors — Seurat SNN Graph Construction
#'
#' PURPOSE:
#'   Constructs a Shared Nearest Neighbor (SNN) graph in PCA (or corrected)
#'   embedding space. Direct prerequisite for FindClusters() and RunUMAP().
#'
#' INPUTS:
#'   input_data  : character (RDS path) | Seurat object
#'   params_list : list (optional, overrides defaults)
#'
#' DEFAULTS:
#'   dims = 1:20, k.param = 20L, reduction = "pca",
#'   prune.SNN = 0.0667, annoy.metric = "euclidean"
#'
#' OUTPUTS:
#'   Seurat object with @graphs$RNA_nn, @graphs$RNA_snn,
#'   and @misc$analysis_history updated

library(Seurat)

run_findneighbors <- function(input_data, params_list = NULL) {
    # Memory-First
    if (is.character(input_data)) {
        obj <- readRDS(input_data)
    } else {
        obj <- input_data
    }

    # Parameter Safety
    default_params <- list(
        dims         = 1:20,
        k.param      = 20L,
        reduction    = "pca",
        prune.SNN    = 0.0667,
        annoy.metric = "euclidean"
    )
    current_params <- modifyList(
        default_params,
        if (!is.null(params_list)) params_list else list()
    )

    # Validate reduction exists
    if (!current_params$reduction %in% names(obj@reductions)) {
        available <- names(obj@reductions)
        stop(paste0(
            "Reduction '", current_params$reduction, "' not found. ",
            "Available: ", paste(available, collapse = ", "), ". ",
            "Run RunPCA() (or integration) first."
        ))
    }

    # Cap dims to available PCs
    max_dim <- ncol(obj@reductions[[current_params$reduction]]@cell.embeddings)
    dims_used <- current_params$dims[current_params$dims <= max_dim]

    obj <- Seurat::FindNeighbors(
        obj,
        dims         = dims_used,
        k.param      = current_params$k.param,
        reduction    = current_params$reduction,
        prune.SNN    = current_params$prune.SNN,
        annoy.metric = current_params$annoy.metric
    )

    n_cells <- ncol(obj)

    # Metadata Footprinting
    if (is.null(obj@misc$analysis_history)) obj@misc$analysis_history <- list()
    obj@misc$analysis_history <- c(obj@misc$analysis_history, list(list(
        skill_id  = "seurat_findneighbors",
        params    = list(
            dims         = as.integer(dims_used),
            k.param      = as.integer(current_params$k.param),
            reduction    = current_params$reduction,
            prune.SNN    = current_params$prune.SNN,
            annoy.metric = current_params$annoy.metric
        ),
        metrics   = list(
            n_dims_used = length(dims_used),
            k_param_used = as.integer(current_params$k.param),
            n_cells     = as.integer(n_cells)
        ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )))

    if (is.character(input_data)) saveRDS(obj, input_data)
    return(obj)
}
