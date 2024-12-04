
#' Extract quality information of features
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' 
#' @return A data frame with quality metrics for each feature.
#'
#' @examples 
#' ex_set <- assess_quality(example_set)
#' quality(ex_set)
#'
#' @export
quality <- function(object) {
  object <- check_object(object)
  if (!all(c("RSD", "RSD_r", "D_ratio", "D_ratio_r") %in%
      colnames(rowData(object)))) {
    return(NULL)
  }
  as.data.frame(
    rowData(object)[c("Feature_ID", "RSD", "RSD_r", "D_ratio", "D_ratio_r")])
}

.erase_quality <- function(object) {
  if (!all(c("RSD", "RSD_r", "D_ratio", "D_ratio_r") %in%
      colnames(rowData(object)))) {
    return(NULL)
  }
  rowData(object)[c("RSD", "RSD_r", "D_ratio", "D_ratio_r")] <- NULL

  object
}



#' Assess quality information of features
#'
#' Assess features using the quality metrics defined in (Broadhurst 
#' 2018). The quality metrics are described in Details section of 
#' \code{flag_quality}
#' @param object a SummarizedExperiment or MetaboSet object
#'
#' @return A SummarizedExperiment or MetaboSet object with quality metrics in 
#' feature data.
#' 
#' @examples
#' ex_set <- assess_quality(example_set)
#' rowData(ex_set)
#'
#' @export
setGeneric("assess_quality", signature = "object",
           function(object, ...) standardGeneric("assess_quality"))

setMethod("assess_quality", signature = c(object = "MetaboSet"),
  function(object) {
    object <- check_object(object, pheno_QC = TRUE, assay.type = 1)
    # Remove old quality metrics
    if (!is.null(quality(object))) {
      object <- .erase_quality(object)
    }

    qc_data <- assay(object)[, object$QC == "QC"]
    sample_data <- assay(object)[, object$QC != "QC"]
    features <- rownames(sample_data)
    
    quality_metrics <- BiocParallel::bplapply(features, function(feature) {
      data.frame(
        Feature_ID = feature,
        RSD = finite_sd(qc_data[feature, ]) / abs(finite_mean(qc_data[feature,])),
        RSD_r = finite_mad(qc_data[feature, ]) /
          abs(finite_median(qc_data[feature, ])),
        D_ratio = finite_sd(qc_data[feature, ]) /
          finite_sd(sample_data[feature, ]),
        D_ratio_r = finite_mad(qc_data[feature, ]) /
          finite_mad(sample_data[feature, ]),
        row.names = feature, stringsAsFactors = FALSE)
    })
  quality_metrics <- do.call(rbind, quality_metrics)
  object <- join_rowData(object, quality_metrics)
})

setMethod("assess_quality", signature = c(object = "SummarizedExperiment"),
  function(object, assay.type = NULL) {
    from <- .get_from_name(object, assay.type)
    object <- check_object(object, pheno_QC = TRUE, assay.type = from)
    # Remove old quality metrics
    if (!is.null(quality(object))) {
      object <- .erase_quality(object)
    }

    qc_data <- assay(object)[, object$QC == "QC"]
    sample_data <- assay(object)[, object$QC != "QC"]
    features <- rownames(sample_data)
    
    quality_metrics <- BiocParallel::bplapply(features, function(feature) {
      data.frame(
        Feature_ID = feature,
        RSD = finite_sd(qc_data[feature, ]) / abs(finite_mean(qc_data[feature,])),
        RSD_r = finite_mad(qc_data[feature, ]) /
          abs(finite_median(qc_data[feature, ])),
        D_ratio = finite_sd(qc_data[feature, ]) /
          finite_sd(sample_data[feature, ]),
        D_ratio_r = finite_mad(qc_data[feature, ]) /
          finite_mad(sample_data[feature, ]),
        row.names = feature, stringsAsFactors = FALSE)
    })
    quality_metrics <- do.call(rbind, quality_metrics)
    object <- join_rowData(object, quality_metrics)
})
  
#' Flag low-quality features
#'
#' Flags low-quality features using the quality metrics defined in (Broadhurst 
#' 2018). The metrics are described in more detain in Details. A condition for 
#' keeping the features is given as a character, which is passed to 
#' \code{dplyr::filter}.
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param condition character, condition for keeping the features, see Details
#'
#' @details The quality metrics measure two things: internal spread of the QCs,
#' and spread of the QCs compared to the spread of the biological samples.
#' Internal spread is measured with relative standard deviation (RSD), also 
#' known as coefficient of variation (CV).
#' \deqn{RSD = sd(QC) / mean(QC) }
#' Where \eqn{sd(QC)} is the standard deviation of the QC samples and 
#' '\eqn{mean(QC)} is the sample mean of the signal in the QC samples.
#' RSD can also be replaced by a non-parametric, robust version based on the 
#' median and median absolute deviation (MAD):
#'   \deqn{RSD_r = 1.4826 * MAD(QC) / median(QC)}
#' The spread of the QC samples compared to the biological samples is measured 
#' using a metric called D-ratio:
#'   \deqn{D_ratio = sd(QC) / sd(biological)}
#' Or, as before, a non-parametric, robust alternative:
#'   \deqn{D_ratio_r = MAD(QC) / MAD(biolofical) }
#' The default condition keeps features that pass either of the two following 
#' conditions:
#' \deqn{RSD_r < 0.2 \& D_ratio_r < 0.4}
#' \deqn{RSD < 0.1 \& RSD_r < 0.1 \& D_ratio < 0.1}
#'
#' @return a SummarizedExperiment or MetaboSet object with the features flagged.
#'
#' @references Broadhurst, David et al. Guidelines and considerations for the 
#' use of system suitability and quality control samples in mass spectrometry 
#' assays applied in untargeted clinical metabolomic studies.
#' Metabolomics : Official journal of the Metabolomic Society vol. 14,6 (2018): 
#' 72. doi:10.1007/s11306-018-1367-3
#'
#' @usage flag_quality
#' flag_quality(object, condition =
#'   "(RSD_r < 0.2 & D_ratio_r < 0.4) | 
#'   (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)")
#'
#' @examples
#' ex_set <- flag_quality(example_set)
#' rowData(ex_set)
#' # Custom condition
#' ex_set <- flag_quality(example_set, 
#'   condition = "RSD_r < 0.3 & D_ratio_r < 0.6")
#' rowData(ex_set)
#'
#' @export
flag_quality <- function(object, assay.type = NULL, 
  condition = "(RSD_r < 0.2 & D_ratio_r < 0.4) | 
              (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)") {
  if (is.null(quality(object))) {
    object <- assess_quality(object, assay.type)
  }
  object <- check_object(object, feature_ID = TRUE)
  .add_citation("Quality metrics were computed as per guidelines in:", paste(
    "Broadhurst, David et al. Guidelines and considerations for the use of",
    "system suitability and quality control samples in mass spectrometry",
    "assays applied in untargeted clinical metabolomic studies.",
    "Metabolomics : Official journal of the Metabolomic Society",
    "vol. 14,6 (2018): 72. doi:10.1007/s11306-018-1367-3"))
  good <- paste0("as.data.frame(rowData(object)) %>% 
                 dplyr::filter(", condition, ")") %>%
    parse(text = .) %>%
    eval()
  good <- good$Feature_ID

  idx <- is.na(flag(object)) & !rowData(object)$Feature_ID %in% good
  flag(object)[idx] <- "Low_quality"

  percentage <- scales::percent(sum(flag(object) == "Low_quality", 
                                    na.rm =TRUE) 
                                / nrow(rowData(object)))
  log_text(paste0("\n", percentage, " of features flagged for low quality"))
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}

#' Flag features with low detection rate
#'
#' Flags features with too high amount of missing values. There are two 
#' detection rate limits, both defined as the minimum proportion of samples 
#' that need to have a value (not NA) for the feature to be kept. 
#' '\code{qc_limit} is the detection rate limit for QC samples, 
#' '\code{group_limit} is the detection rate limit for the actual study groups.
#' If the group limit is passed for AT LEAST ONE GROUP, then the feature is 
#' kept. Features with low detection rate in QCs are flagged as 
#' "Low_qc_detection", while low detection rate in the study groups is flagged 
#' as "Low_group_detection". The detection rates for all the groups are 
#' recorded in feature data.
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param qc_limit the detection rate limit for QC samples
#' @param group_limit the detection rate limit for study groups
#' @param group the columns name in sample information to use as the grouping 
#' variable
#'
#' @return A SummarizedExperiment or MetaboSet object with the features flagged.
#'
#' @examples
#' ex_set <- mark_nas(example_set, value = 0)
#' ex_set <- flag_detection(ex_set, group = "Group")
#' rowData(ex_set)
#'
#' @export
flag_detection <- function(object, qc_limit = 0.7, group_limit = 0.5,
                           group = NULL, assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- check_object(object, pheno_ID = TRUE, pheno_QC = TRUE, 
                         pheno_factors = group, assay.type = from, 
                         feature_ID = TRUE)
  found_qc <- apply(assay(object[, object$QC == "QC"], from), 1, prop_found)
  bad_qc <- names(which(found_qc < qc_limit))

  found_qc_df <- data.frame(Feature_ID = names(found_qc),
                            Detection_rate_QC = found_qc,
                            stringsAsFactors = FALSE)

  idx <- is.na(flag(object)) & rowData(object)$Feature_ID %in% bad_qc
  flag(object)[idx] <- "Low_qc_detection"

  # Compute proportions found in each study group
  if (!is.null(group)) {
    proportions <- combined_data(object)[, c("Sample_ID", group,
                                             rownames(object))] %>%
      tidyr::gather("Feature_ID", "Intensity", rownames(object)) %>%
      dplyr::group_by(.data$Feature_ID, !!as.name(group)) %>%
      dplyr::summarise(proportion_found = prop_found(.data$Intensity)) %>%
      tidyr::spread(!!as.name(group), "proportion_found")
    # Remove a possible QC column
    proportions$QC <- NULL
    colnames(proportions)[-1] <- paste0("Detection_rate_", group, "_",
                                        colnames(proportions)[-1])
    # Check if any group has enough non-missing entries
    proportions$good <- apply(proportions[-1], 1, function(x) {
      any(x >= group_limit)
    })

    idx <- is.na(flag(object)) & (!rowData(object)$Feature_ID %in% 
      proportions$Feature_ID[proportions$good])
    flag(object)[idx] <- "Low_group_detection"
    # Add detection rates to feature data
    proportions <- dplyr::left_join(proportions, found_qc_df, by = "Feature_ID")
    proportions$good <- NULL
  } else {
    proportions <- found_qc_df
  }

  percentage <- scales::percent(sum(flag(object) %in% c("Low_qc_detection",
                                                        "Low_group_detection"), 
                                    na.rm = TRUE) / 
                                nrow(rowData(object)))
  log_text(paste0("\n", percentage, 
                  " of features flagged for low detection rate"))

  object <- join_rowData(object, proportions)
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}


#' Flag contaminants based on blanks
#'
#' Flags contaminant features by comparing the median values of blanks and 
#' biological samples. Biological sampels are defined as samples that are not 
#' marked as blanks and are not QCs. If the median of blanks > the median of 
#' biological samples times a set ratio, the feature is flagged as contaminant.
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param blank_col character, the column name in pheno data with blank labels
#' @param blank_label character, the label for blank samples in blank_col
#' @param flag_thresh numeric, the ratio threshold for flagging contaminants.
#' If the median of blanks > flag_thresh * median of biological samples, the 
#' feature gets flagged.
#' @param flag_label character, the label used when flagging contaminants. Can 
#' be changed if sample processing contaminants and carryover contaminants are 
#' flagged separately.
#'
#' @return A SummarizedExperiment or MetaboSet object with contaminant features 
#' flagged.
#'
#' @examples 
#' # Make a blank sample which has one (first) feature exceeding the threshold
#' ## Abundance matrix
#' med <- median(assay(example_set)[1, example_set$QC != "QC"])
#' assay <- matrix(c(med * 0.05 + 1, rep(0, 79)), ncol = 1, nrow = 80, 
#'                   dimnames = list(NULL, "Demo_51"))
#' assay <- cbind(assay(example_set), assay)
#' ## Sample metadata
#' pheno_data <- colData(example_set)[1, ]
#' rownames(pheno_data) <- "Demo_51"
#' pheno_data$Sample_ID <- "Demo_51"
#' pheno_data$Injection_order <- 51
#' pheno_data[c("Subject_ID", "Group", "QC", "Time")] <- "Blank"
#' pheno_data <- rbind(colData(example_set), pheno_data)
#' ## Feature metadata
#' feature_data <- rowData(example_set)
#'
#' # Construct SummarizedExperiment object with blank sample
#' ex_set <- SummarizedExperiment(assays = assay, 
#'                                colData = pheno_data,
#'                                rowData = feature_data)
#' # Flag contaminant(s)
#' contaminants_flagged <- flag_contaminants(ex_set, blank_col = "QC", 
#'                                           blank_label = "Blank")
#' 
#' @export
flag_contaminants <- function(object, blank_col, blank_label, 
                              flag_thresh = 0.05, flag_label = "Contaminant",
                              assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- check_object(object, pheno_QC = TRUE, pheno_cols = blank_col, 
                         assay.type = from)
  
  blanks <- object[, colData(object)[, blank_col] == blank_label]
  samples <- object[, object$QC != "QC" &
                    colData(object)[, blank_col] != blank_label]

  blank_median <- apply(assay(blanks, from), 1, finite_median)
  sample_median <- apply(assay(samples, from), 1, finite_median)
  blank_flag <- blank_median / sample_median > flag_thresh

  idx <- is.na(flag(object)) & !is.na(blank_flag)
  idx <- idx & blank_flag
  flag(object)[idx] <- flag_label

  percentage <- scales::percent(sum(flag(object) == flag_label, na.rm = TRUE) /
                                nrow(object))
  log_text(paste0("\n", percentage, " of features flagged as contaminants"))

  blank_ratio <- data.frame(Feature_ID = rownames(object), 
                            Blank_ratio = blank_median / sample_median,
                            stringsAsFactors = FALSE)
  object <- join_rowData(object, blank_ratio)
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}
