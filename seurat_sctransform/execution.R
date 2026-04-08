#' seurat_sctransform — SCTransform Normalization
#'
#' PURPOSE:
#'   Regularized negative binomial regression to remove sequencing depth
#'   confounding while preserving biological variance. Replaces the
#'   NormalizeData → FindVariableFeatures → ScaleData pipeline.
#'
#' INPUTS:
#'   input_data  : character (RDS path) | Seurat object
#'   params_list : list (optional, overrides defaults)
#'
#' DEFAULTS:
#'   variable.features.n = 3000, vst.flavor = "v2",
#'   vars.to.regress = NULL, assay = "RNA", verbose = FALSE
#'
#' OUTPUTS:
#'   Seurat object with SCT assay, variable features set,
#'   and @misc$analysis_history updated

library(Seurat)

run_sctransform <- function(input_data, params_list = NULL) {
    # Memory-First
    if (is.character(input_data)) {
        obj <- readRDS(input_data)
    } else {
        obj <- input_data
    }

    # Parameter Safety
    default_params <- list(
        variable.features.n = 3000L,
        vst.flavor          = "v2",
        vars.to.regress     = NULL,
        assay               = "RNA",
        verbose             = FALSE
    )
    current_params <- modifyList(
        default_params,
        if (!is.null(params_list)) params_list else list()
    )

    # Build call args
    sct_args <- list(
        object              = obj,
        variable.features.n = current_params$variable.features.n,
        vst.flavor          = current_params$vst.flavor,
        assay               = current_params$assay,
        verbose             = current_params$verbose
    )
    if (!is.null(current_params$vars.to.regress)) {
        sct_args$vars.to.regress <- current_params$vars.to.regress
    }

    obj <- do.call(Seurat::SCTransform, sct_args)

    # Collect metrics
    n_variable_features <- length(VariableFeatures(obj))
    n_cells <- ncol(obj)

    # Pearson residual variance (median across variable features)
    median_pr_var <- NA_real_
    if ("SCT" %in% names(obj@assays)) {
        sct_model_data <- obj@assays[["SCT"]]@meta.features
        if ("residual_variance" %in% colnames(sct_model_data)) {
            hvf_mask <- rownames(sct_model_data) %in% VariableFeatures(obj)
            median_pr_var <- median(sct_model_data$residual_variance[hvf_mask], na.rm = TRUE)
        }
    }

    # Metadata Footprinting
    if (is.null(obj@misc$analysis_history)) obj@misc$analysis_history <- list()
    obj@misc$analysis_history <- c(obj@misc$analysis_history, list(list(
        skill_id  = "seurat_sctransform",
        params    = list(
            variable.features.n = as.integer(current_params$variable.features.n),
            vst.flavor          = current_params$vst.flavor,
            vars.to.regress     = current_params$vars.to.regress,
            assay               = current_params$assay
        ),
        metrics   = list(
            n_variable_features        = as.integer(n_variable_features),
            median_pearson_residual_var = median_pr_var,
            n_cells                    = as.integer(n_cells)
        ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )))

    if (is.character(input_data)) saveRDS(obj, input_data)
    return(obj)
}
