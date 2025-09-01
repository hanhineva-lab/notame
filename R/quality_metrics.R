
#' Extract quality information of features
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' 
#' @return A data frame with quality metrics for each feature.
#'
#' @examples
#' data(toy_notame_set)
#' ex_set <- assess_quality(toy_notame_set)
#' quality(ex_set)
#'
#' @export
quality <- function(object) {
  object <- .check_object(object)
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
#' \code{\link{flag_quality}}
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A SummarizedExperiment object with quality metrics in 
#' feature data.
#' 
#' @examples
#' data(toy_notame_set)
#' ex_set <- assess_quality(toy_notame_set)
#' rowData(ex_set)
#'
#' @export
assess_quality <- function(object, assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_QC = TRUE, assay.type = from)
  # Remove old quality metrics
  if (!is.null(quality(object))) {
    object <- .erase_quality(object)
  }

  qc_data <- assay(object, from)[, object$QC == "QC"]
  sample_data <- assay(object, from)[, object$QC != "QC"]
  features <- rownames(sample_data)
  
  quality_metrics <- BiocParallel::bplapply(features, function(feature) {
    data.frame(
      Feature_ID = feature,
      RSD = finite_sd(qc_data[feature, ]) / 
      abs(finite_mean(qc_data[feature,])),
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
}

#' Flag low-quality features
#'
#' Flags low-quality features using the quality metrics defined in (Broadhurst 
#' 2018). The metrics are described in more detain in Details. A dual RSD and 
#' D_ratio condition for keeping the features is given.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param RSD character, minimum RSD to keep feature
#' @param D_ratio character, minimum D-ratio to keep feature
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @details The quality metrics measure two things: internal spread of the QCs,
#' and spread of the QCs compared to the spread of the biological samples.
#' Internal spread is measured with relative standard deviation (RSD), also 
#' known as coefficient of variation (CV).
#' RSD used here is the non-parametric, robust version based on the 
#' median and median absolute deviation (MAD):
#'   \deqn{RSD_r = 1.4826 * MAD(QC) / median(QC)}
#' The spread of the QC samples compared to the biological samples is measured 
#' using a metric called D-ratio. The non-parametric, robust alternative is 
#' used here:
#'   \deqn{D_ratio_r = MAD(QC) / MAD(biological) }
#' The default condition keeps features that pass the following condition:
#' \deqn{RSD_r < 0.2 \& D_ratio_r < 0.4}
#'
#' @return a SummarizedExperiment object with the features flagged.
#'
#' @references Broadhurst, David et al. Guidelines and considerations for the 
#' use of system suitability and quality control samples in mass spectrometry 
#' assays applied in untargeted clinical metabolomic studies.
#' Metabolomics : Official journal of the Metabolomic Society vol. 14,6 (2018): 
#' 72. doi:10.1007/s11306-018-1367-3
#'  
#' @examples
#' data(toy_notame_set)
#' ex_set <- flag_quality(toy_notame_set)
#' rowData(ex_set)
#' # Custom condition
#' ex_set <- flag_quality(toy_notame_set, RSD = 0.1, D_ratio = 0.3)
#' rowData(ex_set)
#'
#' @export
flag_quality <- function(object, assay.type = NULL, RSD = 0.2, 
                         D_ratio = 0.4) {
  if (is.null(quality(object))) {
    object <- assess_quality(object, assay.type)
  }
  object <- .check_object(object, feature_ID = TRUE)
  .add_citation("Quality metrics were computed as per guidelines in:", paste(
    "Broadhurst, David et al. Guidelines and considerations for the use of",
    "system suitability and quality control samples in mass spectrometry",
    "assays applied in untargeted clinical metabolomic studies.",
    "Metabolomics : Official journal of the Metabolomic Society",
    "vol. 14,6 (2018): 72. doi:10.1007/s11306-018-1367-3"))

  good_ind <- rowData(object)$RSD_r < RSD & rowData(object)$D_ratio_r < D_ratio
  good <- rowData(object)[good_ind, ]$Feature_ID

  idx <- is.na(flag(object)) & !rowData(object)$Feature_ID %in% good
  flag(object)[idx] <- "Low_quality"

  percentage <- scales::percent(sum(flag(object) == "Low_quality", 
                                    na.rm =TRUE) 
                                / nrow(rowData(object)))
  log_text(paste0("\n", percentage, " of features flagged for low quality"))

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
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param qc_limit the detection rate limit for QC samples
#' @param group_limit the detection rate limit for study groups
#' @param group the columns name in sample information to use as the grouping 
#' variable
#' @param assay.type character, assay to be used in case of multiple assays
#' @return A SummarizedExperiment object with the features flagged.
#'
#' @examples
#' data(toy_notame_set)
#' ex_set <- mark_nas(toy_notame_set, value = 0)
#' ex_set <- flag_detection(ex_set, group = "Group")
#' rowData(ex_set)
#'
#' @export
flag_detection <- function(object, qc_limit = 0.7, group_limit = 0.5,
                           group = NULL, assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_ID = TRUE, pheno_QC = TRUE, 
                         pheno_factors = group, assay.type = from, 
                         feature_ID = TRUE, feature_flag = TRUE)
  found_qc <- apply(assay(object[, object$QC == "QC"], from), 1, prop_found)
  bad_qc <- names(which(found_qc < qc_limit))

  found_qc_df <- data.frame(Feature_ID = names(found_qc),
                            Detection_rate_QC = found_qc,
                            stringsAsFactors = FALSE)

  idx <- is.na(flag(object)) & rowData(object)$Feature_ID %in% bad_qc
  flag(object)[idx] <- "Low_qc_detection"

  # Compute proportions found in each study group
  if (!is.null(group)) {
    proportions <- combined_data(object, from)[, c("Sample_ID", group,
                                               rownames(object))] |>
      tidyr::gather("Feature_ID", "Intensity", rownames(object)) |>
      dplyr::group_by(.data$Feature_ID, !!as.name(group)) |>
      dplyr::summarise(proportion_found = prop_found(.data$Intensity)) |>
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

  object
}


#' Flag contaminants based on blanks
#'
#' Flags contaminant features by comparing the median values of blanks and 
#' biological samples. Biological sampels are defined as samples that are not 
#' marked as blanks and are not QCs. If the median of blanks > the median of 
#' biological samples times a set ratio, the feature is flagged as contaminant.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param blank_col character, the column name in pheno data with blank labels
#' @param blank_label character, the label for blank samples in blank_col
#' @param flag_thresh numeric, the ratio threshold for flagging contaminants.
#' If the median of blanks > flag_thresh * median of biological samples, the 
#' feature gets flagged.
#' @param flag_label character, the label used when flagging contaminants. Can 
#' be changed if sample processing contaminants and carryover contaminants are 
#' flagged separately.
#' @param assay.type character, assay to be used in case of multiple assays
#' @return A SummarizedExperiment object with contaminant features 
#' flagged.
#'
#' @examples
#' data(toy_notame_set)
#' # Make a blank sample which has one (first) feature exceeding the threshold
#' ## Abundance matrix
#' med <- median(assay(toy_notame_set)[1, toy_notame_set$QC != "QC"])
#' assay <- matrix(c(med * 0.05 + 1, rep(0, 79)), ncol = 1, nrow = 80, 
#'                   dimnames = list(NULL, "Demo_51"))
#' assay <- cbind(assay(toy_notame_set), assay)
#' ## Sample metadata
#' pheno_data <- colData(toy_notame_set)[1, ]
#' rownames(pheno_data) <- "Demo_51"
#' pheno_data$Sample_ID <- "Demo_51"
#' pheno_data$Injection_order <- 51
#' pheno_data[c("Subject_ID", "Group", "QC", "Time")] <- "Blank"
#' pheno_data <- rbind(colData(toy_notame_set), pheno_data)
#' ## Feature metadata
#' feature_data <- rowData(toy_notame_set)
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
  object <- .check_object(object, pheno_QC = TRUE, pheno_cols = blank_col, 
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

  object
}
