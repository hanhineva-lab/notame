# Used to combine and return multiple objects from a loop
.comb <- function(x, ...) {
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}

.calc_cubic_spline <- function(feature, sample_data, log_transform, 
                               spar, spar_lower, spar_upper) {
  # Spline cannot be fitted if there are les than 4 QC values
  qc_detected <- sample_data$QC == "QC" & !is.na(feature)
  if (sum(qc_detected) < 4) {
    return(rep(NA, length(feature)))
  }
  # Spline regression on the QC samples
  fit <- stats::smooth.spline(x = sample_data$Injection_order[qc_detected], 
                              y = feature[qc_detected], 
                              all.knots = TRUE, spar = spar, 
                              control.spar = 
                              list("low" = spar_lower, "high" = spar_upper))
  # Predicted values for all samples
  predicted <- stats::predict(fit, sample_data$Injection_order)$y
  # Substraction in log space, division in original space
  if (log_transform) {
    corrected <- feature + 
      mean(feature[qc_detected]) - predicted
  } else {
    corr_factors <- predicted[1] / predicted
    corrected <- feature * corr_factors
  }

  corrected
}

#' Fit a cubic spline to correct drift
#'
#' Corrects the drift in the features by applying smoothed cubic spline 
#' regression to each feature separately.
#'
#' @param object a SummarizedExperiment object
#' @param log_transform logical, should drift correction be done on 
#' log-transformed values? See Details
#' @param spar smoothing parameter
#' @param spar_lower,spar_upper lower and upper limits for the smoothing 
#' parameter
#'
#' @return A list including a SummarizedExperiment object
#' with drift corrected features and predicted = matrix of the predicted values 
#' by the cubic spline (used in visualization).
#'
#' @details If \code{log_transform = TRUE}, the correction will be done on 
#' log-transformed values.
#' The correction formula depends on whether the correction is run on original 
#' values or log-transformed values.
#' In log-space: \eqn{corrected = original + mean of QCs - 
#' prediction by cubic spline}.
#' In original space: \eqn{corrected = original * prediction for first QC / 
#' prediction for current point}.
#' We recommend doing the correction in the log-space since the log-transformed 
#' data better follows the assumptions of cubic spline regression. The drift 
#' correction in the original space also sometimes results in negative values, 
#' and results in rejection of the drift corrrection procedure.
#'
#' If \code{spar} is set to \code{NULL} (the default), the smoothing parameter 
#' will be separately chosen for each feature from the range 
#' [\code{spar_lower, spar_upper}] using cross validation.
#'
#' @examples
#' data(toy_notame_set)
#' dc <- dc_cubic_spline(toy_notame_set)
#' corrected <- dc$object
#'
#' @seealso  \code{\link[stats]{smooth.spline}} for details about the 
#' regression, \code{\link{inspect_dc}} for analysing the drift correction 
#' results, \code{\link{save_dc_plots}} for plotting the drift correction 
#' process for each feature
#'
#' @noRd
dc_cubic_spline <- function(object, log_transform = TRUE, spar = NULL,
                            spar_lower = 0.5, spar_upper = 1.5, 
                            assay.type = NULL, name = NULL) {
  # Start log
  log_text("Starting drift correction")
  # Zero values do not behave correctly
  full_data <- assay(object, assay.type)
  if (sum(full_data == 0, na.rm = TRUE)) {
    log_text(paste0("Zero values in feature abundances detected.",
                    " Zeroes will be replaced with 1.1."))
    full_data[full_data == 0] <- 1.1
  }
  
  # log-transform before fiting the cubic spline
  if (log_transform) {
    full_data <- log(full_data)
  }

  # Return both predicted values (for plotting) and drift corrected values
  dc_data <- BiocParallel::bplapply(
    as.data.frame(t(full_data)),
    .calc_cubic_spline,
    colData(object),
    log_transform,
    spar,
    spar_lower,
    spar_upper
  )
  corrected <- do.call(rbind, dc_data)
  colnames(corrected) <- colnames(full_data)
  # Inverse the initial log transformation
  if (log_transform) {
    corrected <- exp(corrected)
  }
  assay(object, name) <- corrected
  # Recompute quality metrics
  log_text("Recomputing quality metrics for drift corrected data")
  object <- assess_quality(object, assay.type = name)
  
  log_text("Drift correction performed")

  object
}


.help_inspect_dc <- function(dc, orig, check_quality, condition,
                             assay.orig, assay.dc) {
  orig_data <- assay(orig, assay.orig)
  dc_data <- assay(dc, assay.dc)
  features <- rownames(orig)
  rowData(dc)$DC_note <- NA
  missing_qcs <- apply(dc_data, 1, \(x) all(is.na(x)))
  rowData(dc)$DC_note[missing_qcs] <- "Missing_QCS"
  negative_dc <- apply(dc_data, 1, \(x) any(x < 0, na.rm = TRUE))
  rowData(dc)$DC_note[negative_dc & !missing_qcs] <- "Negative_DC"
  not_passing <- rep(FALSE, length(features))
  if (check_quality){
    qdiff <- quality(dc)[2:5] - quality(orig)[2:5] |> 
      as.data.frame() |> 
      t()
    pass <- paste0("qdiff |> dplyr::filter(", condition, ")") |>
      parse(text = _) |>
      eval() |>
      rownames()
    not_passing <- (features %in% !pass) & !(missing_qcs | negative_dc)
    rowData(dc)$DC_note[not_passing] <- "Low_quality"
  }
  dc_passing <- !(missing_qcs | negative_dc | not_passing)
  rowData(dc)$DC_note[dc_passing] <- "Drift_corrected"
  #  Replace the features that did not pass with the original values
  assay(dc, assay.dc)[!dc_passing, ] <- orig_data[!dc_passing, ]
  quality_cols <- c("RSD", "RSD_r", "D_ratio", "D_ratio_r")
  rowData(dc)[!dc_passing, quality_cols] <-
    rowData(orig)[!dc_passing, quality_cols]
  dc
}

#' Flag the results of drift correction
#'
#' Determines whether the drift correction worked.
#' The primary reason is to search for features where there were too many 
#' missing values in the QCs, so it was not possible to run drift correction. 
#' If the drift correction is run on the original values (not log-transformed), 
#' then there is also a need to check that the correction did not result
#' in any negative values. This can sometimes happen if the prediction curve 
#' takes an extreme shape.
#'
#' If quality is monitored, a quality condition is checked for each feature. 
#' If the condition is fulfilled, the drift corrected feature is retained,
#' otherwise the original feature is retained and the drift corrected feature 
#' is discarded. The result of this operation is recorded in the feature data.
#'
#' @param orig a SummarizedExperiment object, before drift correction
#' @param dc a SummarizedExperiment object, after drift correction
#' @param check_quality logical, whether quality should be monitored.
#' @param condition a character specifying the condition, see Details
#'
#' @return A SummarizedExperiment object.
#'
#' @details The \code{condition} parameter should be a character giving a 
#' condition compatible with \code{\link[dplyr]{filter}}. The condition is 
#' applied on the \strong{changes} in the quality metrics RSD, RSD_r, D_ratio 
#' and D_ratio_r. 
#' For example, the default is "RSD_r < 0 and D_ratio_r < 0",
#' meaning that both RSD_r and D_ratio_r need to decrease in the drift 
#' correction, otherwise the  drift corrected feature is discarded and the 
#' original is retained.
#'
#' @seealso \code{\link{correct_drift}}, \code{\link{save_dc_plots}}
#'
#' @examples
#' data(toy_notame_set)
#' dc <- dc_cubic_spline(toy_notame_set)
#' corrected <- dc$object
#' inspected <- inspect_dc(
#'   orig = toy_notame_set, dc = corrected,
#'   check_quality = TRUE
#' )
#'
#' @noRd
inspect_dc <- function(orig, dc, check_quality,
                       condition = "RSD_r < 0 & D_ratio_r < 0",
                       assay.orig = NULL, assay.dc = NULL) {
  log_text("Inspecting drift correction results")
  if (is.null(quality(orig))) {
    log_text("Original quality metrics missing, recomputing")
    orig <- assess_quality(orig, assay.type = assay.orig)
  }
  if (is.null(quality(dc))) {
    log_text("Drift corrected quality metrics missing, recomputing")
    dc <- assess_quality(dc, assay.type = assay.dc)
  }
  
  dc <- .help_inspect_dc(dc, orig, check_quality, condition,
                         assay.orig, assay.dc)

  # Log information
  dc_note <- rowData(dc)$DC_note
  note_counts <- table(dc_note) |> unname()
  note_percentage <- note_counts / sum(note_counts)
  note_percentage <- scales::percent(as.numeric(note_percentage))
  note_labels <- table(dc_note) |> names()
  report <- paste(note_labels, note_percentage, sep = ": ", collapse = ",  ")
  log_text(paste0("Drift correction results inspected: ", report))

  dc
}

#' Correct drift using cubic spline
#'
#' A wrapper function for applying cubic spline drift correction and saving
#' before and after plots.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param log_transform logical, should drift correction be done on 
#' log-transformed values? See Details
#' @param spar smoothing parameter as in 
#' \code{\link[stats]{smooth.spline}}
#' @param spar_lower,spar_upper lower and upper limits for the smoothing 
#' parameter
#' @param check_quality logical, whether quality should be monitored.
#' @param condition a character specifying the condition used to decide whether 
#' drift correction works adequately, see Details
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param color_scale,shape_scale the color and shape scales as returned by a 
#' ggplot function
#' @param assay.type character, assay to be used in case of multiple assays
#' @param name character, name of the resultant assay
#'
#' @return A SummarizedExperiment object as the one supplied, with 
#' drift corrected features.
#'
#' @details If \code{log_transform = TRUE}, the correction will be done on 
#' log-transformed values.
#' The correction formula depends on whether the correction is run on original 
#' values or log-transformed values.
#' In log-space: \eqn{corrected = original + mean of QCs - prediction by cubic 
#' spline}.
#' In original space: \eqn{corrected = original * prediction for first QC / 
#' prediction for current point}.
#' We recommend doing the correction in the log-space since the log-transfomred 
#' data better follows the assumptions of cubic spline regression. The drift 
#' correction in the original space also sometimes results
#' in negative values, and results in rejection of the drift corrrection 
#' procedure.
#' If  \code{check_quality = TRUE}, the \code{condition} parameter should be a 
#' character giving a condition compatible with \code{\link[dplyr]{filter}}. 
#' The condition is applied on the \strong{changes} in the quality metrics
#' RSD, RSD_r, D_ratio and D_ratio_r. For example, the default is "RSD_r < 0 
#' and D_ratio_r < 0", meaning that both RSD_r and D_ratio_r need to decrease
#' in the drift correction, otherwise the drift corrected feature is discarded 
#' and the original is retained.
#' By default, the column used for color is also used for shape.
#'
#' @examples
#' data(toy_notame_set)
#' corrected <- correct_drift(mark_nas(toy_notame_set[1:5, ], value = 0))
#'
#' @seealso \code{\link[stats]{smooth.spline}} for details about the regression
#'
#' @export
correct_drift <- function(object, log_transform = TRUE, spar = NULL, 
                          spar_lower = 0.5, spar_upper = 1.5,
                          check_quality = FALSE, 
                          condition = "RSD_r < 0 & D_ratio_r < 0",
                          file = NULL, width = 16, 
                          height = 8, color = "QC", shape = color, 
                          color_scale = getOption("notame.color_scale_dis"),
                          shape_scale = scale_shape_manual(values = c(15, 16)),
                          assay.type = NULL, name = NULL) {
  from_to <- .get_from_to_names(object, assay.type, name)
  object <- .check_object(object, pheno_injection = TRUE, pheno_QC = TRUE, 
                          assay.type = from_to[[1]])
                          
  # Fit cubic spline and correct
  corrected <- dc_cubic_spline(object, log_transform = log_transform, 
                               spar = spar, spar_lower = spar_lower,
                               spar_upper = spar_upper, 
                               assay.type = from_to[[1]],
                               name = from_to[[2]])
  # Only keep corrected versions of features with increased quality
  inspected <- inspect_dc(orig = object, dc = corrected, 
                          check_quality = check_quality, condition = condition, 
                          assay.orig = from_to[[1]], assay.dc = from_to[[2]])
  # Return the final version
  inspected
}