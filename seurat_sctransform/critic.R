#' seurat_sctransform critic — Evaluates SCTransform normalization quality
#'
#' FEEDBACK KEYS (match parameter_science_guide in skill.yaml):
#'   too_few_hvf          → variable.features.n too low
#'   noisy_features       → variable.features.n too high
#'   cell_cycle_confounded → cell cycle covariates not regressed
#'   mito_confounded      → mitochondrial content not regressed

critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol             <- if (!is.null(context$protocol)) context$protocol else "10x Genomics"
    regress_cell_cycle   <- if (!is.null(context$regress_cell_cycle)) context$regress_cell_cycle else FALSE

    thresholds <- list(
        "10x Genomics" = list(min_hvf = 1000L, max_hvf = 5000L, min_pr_var = 0.5),
        "Smart-seq2"   = list(min_hvf = 2000L, max_hvf = 8000L, min_pr_var = 0.4),
        "Drop-seq"     = list(min_hvf = 1000L, max_hvf = 4000L, min_pr_var = 0.5),
        "CITE-seq"     = list(min_hvf = 1500L, max_hvf = 5000L, min_pr_var = 0.4)
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]

    metrics  <- list()
    warnings <- character(0)
    feedback <- character(0)

    # Check SCT assay exists
    if (!"SCT" %in% names(obj@assays)) {
        return(list(
            metrics      = list(),
            warnings     = "CRITICAL: SCT assay not found — run SCTransform() first",
            feedback     = character(0),
            success      = FALSE,
            context_used = protocol
        ))
    }

    n_hvf <- length(VariableFeatures(obj))
    metrics$n_variable_features <- as.integer(n_hvf)
    metrics$n_cells <- as.integer(ncol(obj))

    # Pearson residual variance
    sct_meta <- obj@assays[["SCT"]]@meta.features
    if ("residual_variance" %in% colnames(sct_meta)) {
        hvf_mask <- rownames(sct_meta) %in% VariableFeatures(obj)
        metrics$median_pearson_residual_var <- median(sct_meta$residual_variance[hvf_mask], na.rm = TRUE)
    }

    if (n_hvf < t$min_hvf) {
        warnings <- c(warnings, sprintf(
            "TOO FEW VARIABLE FEATURES: %d (min %d for %s). Increase variable.features.n.",
            n_hvf, t$min_hvf, protocol
        ))
        feedback <- c(feedback, "too_few_hvf")
    }

    if (n_hvf > t$max_hvf) {
        warnings <- c(warnings, sprintf(
            "TOO MANY VARIABLE FEATURES: %d (max %d for %s). Consider decreasing variable.features.n.",
            n_hvf, t$max_hvf, protocol
        ))
        feedback <- c(feedback, "noisy_features")
    }

    # Check if cell cycle scores exist but not regressed
    has_s_score  <- "S.Score"   %in% colnames(obj@meta.data)
    has_g2m      <- "G2M.Score" %in% colnames(obj@meta.data)
    if (regress_cell_cycle && (has_s_score || has_g2m)) {
        # Check if they were regressed
        history <- obj@misc$analysis_history
        regressed_cc <- FALSE
        if (!is.null(history)) {
            for (entry in rev(history)) {
                if (!is.null(entry$skill_id) && entry$skill_id == "seurat_sctransform") {
                    vtr <- entry$params$vars.to.regress
                    if (!is.null(vtr) && any(c("S.Score", "G2M.Score") %in% vtr)) {
                        regressed_cc <- TRUE
                    }
                    break
                }
            }
        }
        if (!regressed_cc) {
            warnings <- c(warnings,
                "CELL CYCLE NOT REGRESSED: S.Score/G2M.Score present but not in vars.to.regress. ",
                "Cell cycle may drive clustering. Add vars.to.regress = c('S.Score', 'G2M.Score')."
            )
            feedback <- c(feedback, "cell_cycle_confounded")
        }
    }

    # Check mitochondrial content
    mito_col <- intersect(c("percent.mt", "pct_counts_mt", "percent_mito"), colnames(obj@meta.data))
    if (length(mito_col) > 0) {
        mito_vals <- obj@meta.data[[mito_col[1]]]
        if (median(mito_vals, na.rm = TRUE) > 15) {
            warnings <- c(warnings, sprintf(
                "HIGH MITOCHONDRIAL CONTENT: median %.1f%%. Consider adding 'percent.mt' to vars.to.regress.",
                median(mito_vals, na.rm = TRUE)
            ))
            feedback <- c(feedback, "mito_confounded")
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
