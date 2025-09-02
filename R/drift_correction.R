# Used to combine and return multiple objects from a loop
.comb <- function(x, ...) {
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}

.calc_cubic_spline <- function(feature, full_data, full_order, qc_data, n,
                               qc_order, log_transform, 
                               spar, spar_lower, spar_upper) {
  dnames <- list(feature, colnames(full_data))
  # Spline cannot be fitted if there are les than 4 QC values
  qc_detected <- !is.na(qc_data[feature, ])
  if (sum(qc_detected) < 4) {
    return(list(
      corrected = matrix(NA_real_, nrow = 1, ncol = n, dimnames = dnames),
      predicted = matrix(NA_real_, nrow = 1, ncol = n, dimnames = dnames)))
  }
  # Spline regression
  fit <- stats::smooth.spline(x = qc_order[qc_detected], 
                              y = qc_data[feature, qc_detected], 
                              all.knots = TRUE, spar = spar, 
                              control.spar = 
                              list("low" = spar_lower, "high" = spar_upper))
  predicted <- stats::predict(fit, full_order)$y
  # Substraction in log space, division in original space
  if (log_transform) {
    corrected <- full_data[feature, ] + 
      mean(qc_data[feature, qc_detected]) - predicted
  } else {
    corr_factors <- predicted[1] / predicted
    corrected <- full_data[feature, ] * corr_factors
  }
  # Each iteration of the loop returns one row to both corrected and predicted
  list(corrected = matrix(corrected, ncol = n, dimnames = dnames),
       predicted = matrix(predicted, ncol = n, dimnames = dnames))
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
                            assay.type = NULL, name = NULL,
                            name_predicted = NULL) {
  # Start log
  log_text(paste("\nStarting drift correction at", Sys.time()))
  # Zero values do not behave correctly
  full_data <- assay(object, assay.type)
  if (sum(full_data == 0, na.rm = TRUE)) {
    log_text(paste0("Zero values in feature abundances detected.",
                    " Zeroes will be replaced with 1.1."))
    full_data[full_data == 0] <- 1.1
  }
  # Extract data and injection order for QC samples and the full dataset
  features <- rownames(object)
  qc_data <- full_data[, object$QC == "QC"]
  qc_order <- object[, object$QC == "QC"]$Injection_order
  full_order <- object$Injection_order
  
  # log-transform before fiting the cubic spline
  if (log_transform) {
    qc_data <- log(qc_data)
    full_data <- log(full_data)
  }
  # comb needs a matrix with the right amount of columns
  n <- ncol(full_data)
  # Return both predicted values (for plotting) and drift corrected values
  dc_data <- BiocParallel::bplapply(features, .calc_cubic_spline, full_data,
                                    full_order, qc_data, n, qc_order,
                                    log_transform, spar, spar_lower, spar_upper)
  dc_data <- do.call(.comb, dc_data)
  corrected <- dc_data$corrected
  # Inverse the initial log transformation
  if (log_transform) {
    corrected <- exp(corrected)
  }
  assay(object, name_predicted) <- dc_data$predicted
  assay(object, name) <- corrected
  # Recompute quality metrics
  object <- assess_quality(object, assay.type = name)
  
  log_text(paste("Drift correction performed at", Sys.time()))

  object
}


.help_inspect_dc <- function(feature, orig_data, dc_data, qdiff,
                             check_quality, condition) {
  data <- orig_data[feature, ]
  if (all(is.na(dc_data[feature, ]))) {
    dc_note <- "Missing_QCS"
  } else if (any(dc_data[feature, ] < 0, na.rm = TRUE)) {
    dc_note <- "Negative_DC"
  } else if (check_quality) {
    pass <- paste0("qdiff[feature, ] |> dplyr::filter(", condition, 
                   ") |> nrow() |> as.logical()") |>
      parse(text = _) |>
      eval()
    if (!pass) {
      dc_note <- "Low_quality"
    } else {
      data <- dc_data[feature, ]
      dc_note <- "Drift_corrected"
    }
  } else {
    data <- dc_data[feature, ]
    dc_note <- "Drift_corrected"
  }

  list(data = matrix(data, nrow = 1, dimnames = list(feature, names(data))),
       dc_notes = data.frame(Feature_ID = feature, DC_note = dc_note,
                             stringsAsFactors = FALSE))
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
                       assay.orig = NULL, assay.dc = NULL, name = NULL) {
  if (is.null(quality(orig))) {
    orig <- assess_quality(orig, assay.type = assay.orig)
  }
  if (is.null(quality(dc))) {
    dc <- assess_quality(dc, assay.type = assay.dc)
  }

  orig_data <- assay(orig, assay.orig)
  dc_data <- assay(dc, assay.dc)
  features <- rownames(orig)
  qdiff <- quality(dc)[2:5] - quality(orig)[2:5]

  log_text(paste("Inspecting drift correction results", Sys.time()))

  inspected <- BiocParallel::bplapply(features, .help_inspect_dc, orig_data,
                                      dc_data, qdiff, check_quality, condition)
  
  inspected <- do.call(.comb, inspected)

  assay(dc, name) <- inspected$data
  dc <- assess_quality(dc, assay.type = name)
  dc <- join_rowData(dc, inspected$dc_notes)

  log_text(paste("Drift correction results inspected at", Sys.time()))

  # Log information
  dc_note <- inspected$dc_notes$DC_note
  note_counts <- table(dc_note) |> unname()
  note_percentage <- note_counts / sum(note_counts)
  note_percentage <- scales::percent(as.numeric(note_percentage))
  note_labels <- table(dc_note) |> names()
  report <- paste(note_labels, note_percentage, sep = ": ", collapse = ",  ")
  log_text(paste0("\nDrift correction results inspected, report:\n", report))

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
#' @param RSD character, change in RSD to keep feature
#' @param D_ratio character, change in D_ratio to keep feature
#' @param plotting logical, whether plots should be drawn
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param color_scale,shape_scale the color and shape scales as returned by a 
#' ggplot function
#' @param assay.type character, assay to be used in case of multiple assays
#' @param name character, name of the resultant assay
#' @param name_predicted character, name of the resultant assay with predicted 
#' values (for plotting)
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
#' If \code{spar} is set to \code{NULL} (the default), the smoothing parameter 
#' will be separately chosen for each feature from the range
#' [\code{spar_lower, spar_upper}] using cross validation. 
#' If  \code{check_quality = TRUE}, the RSD and D_ratio condition is applied on 
#' the \strong{changes} in the quality metrics RSD_r, and D_ratio_r. For 
#' example, the default is RSD_r < 0 and D_ratio_r < 0, meaning that both RSD_r 
#' and D_ratio_r need to decrease 
#' in the drift correction, otherwise the drift corrected feature is discarded 
#' and the original is retained.
#' By default, the column used for color is also used for shape.
#'
#' @examples
#' data(toy_notame_set)
#' corrected <- correct_drift(mark_nas(toy_notame_set[1:5, ], value = 0),
#'   file = "drift_plots.pdf", plotting = TRUE)
#'
#' @seealso \code{\link[stats]{smooth.spline}} for details about the regression
#'
#' @export
correct_drift <- function(object, log_transform = TRUE, spar = NULL, 
                          spar_lower = 0.5, spar_upper = 1.5,
                          check_quality = FALSE, 
                          RSD = 0, D_ratio = 0, 
                          plotting = FALSE, file = NULL, width = 16, 
                          height = 8, color = "QC", shape = color, 
                          color_scale = getOption("notame.color_scale_dis"),
                          shape_scale = scale_shape_manual(values = c(15, 16)),
                          assay.type = NULL, name = "corrected",
                          name_predicted = "drift_pred") {
  from_to <- .get_from_to_names(object, assay.type, name)
  object <- .check_object(object, pheno_injection = TRUE, pheno_QC = TRUE, 
                          assay.type = from_to[[1]])
                          
  # Fit cubic spline and correct
  corrected <- dc_cubic_spline(object, log_transform = log_transform, 
                               spar = spar, spar_lower = spar_lower,
                               spar_upper = spar_upper, 
                               assay.type = from_to[[1]],
                               name = from_to[[2]],
                               name_predicted = name_predicted)
  # Only keep corrected versions of features with increased quality
  condition = paste0("RSD_r < ", RSD, " & D_ratio_r < ", D_ratio)
  inspected <- inspect_dc(orig = object, dc = corrected, 
                          check_quality = check_quality, condition = condition, 
                          assay.orig = from_to[[1]],
                          assay.dc = from_to[[2]], name = from_to[[2]])
  # Optionally save before and after plots
  if (plotting) {
    if (is.null(file)) {
      stop("File must be specified.")
    }
    if (!requireNamespace("notameViz", quietly = TRUE)) {
      stop("Package \"notameViz\" needed for this function to work.",
           " Please install it.", call. = FALSE)
    }
    notameViz::save_dc_plots(inspected, file = file,
      log_transform = log_transform, width = width, height = height,
      color = color, shape = shape, color_scale = color_scale,
      shape_scale = shape_scale, assay.orig = from_to[[1]],
      assay.dc = from_to[[2]], assay.pred = name_predicted)
  }
  # Return the final version
  inspected
}