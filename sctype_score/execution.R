#' sctype_score — ScType Marker Gene Scoring for Cell Type Identification
#'
#' PURPOSE:
#'   Automatically assigns cell types to clusters using marker gene enrichment
#'   scoring. Uses either built-in tissue-specific marker databases or custom
#'   XLSX files. Requires scaled expression data.
#'
#' INPUTS:
#'   input_data  : character (RDS path) | Seurat object
#'   params_list : list (optional, overrides defaults)
#'
#' DEFAULTS:
#'   tissue_type = "Immune system", db_path = NULL, scaled_assay = "scale.data",
#'   cluster_col = "seurat_clusters", score_col = "sctype_classification"
#'
#' OUTPUTS:
#'   Seurat object with ScType classifications in @meta.data
#'   and @misc$analysis_history updated

# ScType source functions (loaded from GitHub or local cache)
.load_sctype <- function() {
    # Try to load from installed package first
    if (requireNamespace("sctype", quietly = TRUE)) {
        return(invisible(NULL))
    }
    # Load core functions from GitHub
    sctype_score_url  <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R"
    gene_sets_url     <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"
    for (url in c(sctype_score_url, gene_sets_url)) {
        source(url)
    }
}

run_sctype_score <- function(input_data, params_list = NULL) {
    # Memory-First
    if (is.character(input_data)) {
        obj <- readRDS(input_data)
    } else {
        obj <- input_data
    }

    # Parameter Safety
    default_params <- list(
        tissue_type    = "Immune system",
        db_path        = NULL,
        scaled_assay   = "scale.data",
        cluster_col    = "seurat_clusters",
        score_col      = "sctype_classification",
        top_cell_types = 1L
    )
    current_params <- modifyList(
        default_params,
        if (!is.null(params_list)) params_list else list()
    )

    # Load ScType functions
    .load_sctype()

    # Resolve marker database
    db <- if (!is.null(current_params$db_path)) {
        current_params$db_path
    } else {
        # Built-in ScType database from GitHub
        "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    }

    # Prepare gene sets
    gs_list <- gene_sets_prepare(db, current_params$tissue_type)

    # Get scaled expression matrix
    assay_name <- current_params$scaled_assay
    if (assay_name == "scale.data") {
        scaled_mat <- Seurat::GetAssayData(obj, slot = "scale.data")
    } else if (assay_name == "SCT") {
        scaled_mat <- Seurat::GetAssayData(obj, assay = "SCT", slot = "scale.data")
    } else {
        scaled_mat <- Seurat::GetAssayData(obj, assay = assay_name, slot = "scale.data")
    }

    if (nrow(scaled_mat) == 0) {
        stop(paste0(
            "Scaled data not found in assay '", assay_name, "'. ",
            "Run ScaleData() or SCTransform() first."
        ))
    }

    # Compute ScType scores
    es_max <- sctype_score(
        scRNAseqData = scaled_mat,
        scaled       = TRUE,
        gs           = gs_list$gs_positive,
        gs2          = gs_list$gs_negative
    )

    # Assign cell types per cluster
    cluster_col <- current_params$cluster_col
    if (!cluster_col %in% colnames(obj@meta.data)) {
        # Fallback to available cluster column
        for (alt in c("seurat_clusters", "leiden", "louvain")) {
            if (alt %in% colnames(obj@meta.data)) {
                cluster_col <- alt
                break
            }
        }
    }

    clusters   <- obj@meta.data[[cluster_col]]
    cL_results <- do.call("rbind", lapply(unique(clusters), function(cl) {
        es_cl <- sort(rowSums(es_max[, rownames(obj@meta.data)[clusters == cl], drop = FALSE]),
                      decreasing = TRUE)
        n_cells_cl <- sum(clusters == cl)
        head(data.frame(
            cluster   = cl,
            type      = names(es_cl),
            scores    = es_cl / n_cells_cl,
            ncells    = n_cells_cl
        ), current_params$top_cell_types)
    }))

    # Best type per cluster
    best_types <- cL_results %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(1, scores) %>%
        dplyr::mutate(
            type = ifelse(scores < 0, "Unknown", as.character(type))
        )

    # Map back to cells
    cell_types <- setNames(
        as.character(best_types$type[match(clusters, best_types$cluster)]),
        rownames(obj@meta.data)
    )

    score_col <- current_params$score_col
    obj@meta.data[[score_col]] <- cell_types

    # Per-cell cluster score
    cell_cluster_scores <- cL_results$scores[match(clusters, cL_results$cluster)]
    obj@meta.data[[paste0(score_col, ".score")]] <- cell_cluster_scores

    # Collect metrics
    n_assigned  <- length(unique(cell_types[cell_types != "Unknown"]))
    pct_unknown <- mean(cell_types == "Unknown") * 100
    mean_score  <- mean(cell_cluster_scores[cell_types != "Unknown"], na.rm = TRUE)
    top3        <- names(sort(table(cell_types[cell_types != "Unknown"]), decreasing = TRUE))[1:min(3L, n_assigned)]

    # Metadata Footprinting
    if (is.null(obj@misc$analysis_history)) obj@misc$analysis_history <- list()
    obj@misc$analysis_history <- c(obj@misc$analysis_history, list(list(
        skill_id  = "sctype_score",
        params    = list(
            tissue_type  = current_params$tissue_type,
            scaled_assay = current_params$scaled_assay,
            cluster_col  = cluster_col
        ),
        metrics   = list(
            n_assigned_types = as.integer(n_assigned),
            pct_unknown      = pct_unknown,
            mean_top_score   = mean_score,
            top3_types       = top3,
            n_clusters       = as.integer(length(unique(clusters))),
            n_cells          = as.integer(ncol(obj))
        ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )))

    if (is.character(input_data)) saveRDS(obj, input_data)
    return(obj)
}
