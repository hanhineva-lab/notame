#' Mark specified values as missing
#'
#' Replaces all values in the peak table that equal the specified value 
#' with NA.
#' For example, vendor software might use 0 or 1 to signal a missing value,
#' which is not understood by R.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param value the value to be converted to NA
#' @param assay.type character, assay to be used in case of multiple assays
#' @param name character, name of the resultant assay in case of multiple assays
#'
#' @return SummarizedExperiment or MetaboSet object as the one supplied, with 
#' missing values correctly set to NA.
#'
#' @examples
#' data(example_set)
#' nas_marked <- mark_nas(example_set, value = 0)
#'
#' @export
mark_nas <- function(object, value, assay.type = NULL, name = NULL) {
  from_to <- .get_from_to_names(object, assay.type, name)
  object <- .check_object(object, assay.type = from_to[[1]])
  ex <- assay(object, from_to[[1]])
  ex[ex == value] <- NA
  assay(object, from_to[[2]]) <- ex
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}

#' Transform the MS/MS output to publication ready
#'
#' Change the MS/MS output from MS-DIAL format to publication-ready format.
#' Original spectra is sorted according to abundance percentage and clarified. 
#' See the example below.
#'
#' Original MS/MS spectra from MS-DIAL:
#' m/z:Raw Abundance
#'
#' 23.193:254 26.13899:5 27.50986:25 55.01603:82 70.1914:16 73.03017:941 
#' 73.07685:13 73.13951:120
#'
#' Spectra after transformation:
#' m/z  (Abundance)
#'
#' 73.03 (100), 23.193 (27), 73.14 (12.8), 55.016 (8.7)
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param ms_ms_spectrum_col name of column with original MS/MS spectra
#' @param peak_num maximum number of peak that is kept (Recommended: 4-10)
#' @param min_abund minimum relative abundance to be kept (Recommended: 1-5)
#' @param deci_num maximum number of decimals to m/z value (Recommended: >2)
#'
#' @return A SummarizedExperiment or MetaboSet object as the one supplied, with 
#' publication-ready MS/MS peak information.
#'
#' @examples
#' data(example_set)
#' # Spectra before fixing
#' ex_set <- example_set
#' rowData(ex_set)$MS_MS_spectrum <- NA
#' rowData(ex_set)[1, ]$MS_MS_spectrum <- 
#'   "28.769:53 44.933:42 52.106:89 69.518:140"
#' rowData(ex_set)$MS_MS_spectrum[
#'   !is.na(rowData(ex_set)$MS_MS_spectrum)]
#' # Fixing spectra with default settings
#' fixed_MSMS_peaks <- fix_MSMS(ex_set)
#' # Spectra after fixing
#' rowData(fixed_MSMS_peaks)$MS_MS_Spectrum_clean[
#'   !is.na(rowData(fixed_MSMS_peaks)$MS_MS_Spectrum_clean)]
#'
#' @export
fix_MSMS <- function(object, ms_ms_spectrum_col = "MS_MS_spectrum",
                     peak_num = 10, min_abund = 5, deci_num = 3) {
  object <- .check_object(object, feature_cols = ms_ms_spectrum_col)
  spec <- rowData(object)[, ms_ms_spectrum_col]
  to_metab <- NULL

  for (i in seq_along(spec)) {
    # Check the feature spectra and skip if it doesn't exist
    if (is.na(spec[i])) {
      to_metab[i] <- NA
      next()
    }
    # Transform format
    spectrum <- spec[i]
    spectrum2 <- t(stringr::str_split(spectrum, pattern = " ", simplify = TRUE))
    spectrum3 <- as.data.frame(stringr::str_split(spectrum2, 
                                                  pattern = ":", 
                                                  simplify = TRUE))
    spectrum3 <- as.data.frame(lapply(spectrum3, as.numeric))
    spectrum3 <- spectrum3[order(spectrum3$V2, decreasing = TRUE), ]

    # Leave n most intense fragment peaks or all peaks if number of peaks < n
    ifelse(nrow(spectrum3) > peak_num, num <- peak_num, num <- nrow(spectrum3))
    spectrum4 <- spectrum3[c(seq_len(num)), ]

    # Round m/z of fragments to n decimals and calculate relative intensity (%)
    spectrum4$V1 <- round(spectrum4$V1, digits = deci_num)
    spectrum4$relative <- round(spectrum4$V2 / max(spectrum3$V2) * 100, 
                                digits = 1)

    # Remove peaks with relative intensity less than n% (recommended: 1-5)
    spectrum5 <- spectrum4[c(spectrum4$relative > min_abund), ]

    # Finalize format and write results
    to_metab[i] <- paste(paste0(spectrum5$V1, " (", spectrum5$relative, ")"),
                         collapse = ", ")
  }

  rowData(object)$MS_MS_Spectrum_clean <- to_metab
  log_text(paste0("Saving fixed MS/MS spectra to column",
                  " \'MS_MS_Spectrum_clean\' in rowData"))
                  
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}

#' Drop QC samples
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#'
#' @return A SummarizedExperiment or MetaboSet object as the one supplied, 
#' without QC samples.
#'
#' @examples
#' data(example_set)
#' dim(example_set)
#' noqc <- drop_qcs(example_set)
#' dim(noqc)
#'
#' @export
drop_qcs <- function(object) {
  object <- .check_object(object)
  object <- object[, object$QC != "QC"]
  colData(object) <- droplevels(colData(object))
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}


#' Drop flagged features
#'
#' Removes all features that have been flagged by quality control functions.
#' Only features that do not have a flag (Flag == NA) are retained.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param all_features logical, should all features be retained? Mainly used by 
#' internal functions
#' 
#' @return A SummarizedExperiment or MetaboSet object without the previously 
#' flagged features.
#'
#' @examples
#' data(example_set)
#' dim(example_set)
#' flagged <- flag_quality(example_set)
#' noflags <- drop_flagged(flagged)
#' dim(noflags)
#'
#' @export
setGeneric("drop_flagged", signature = "object",
           function(object, all_features = FALSE) 
           standardGeneric("drop_flagged"))

#' @rdname drop_flagged
#' @export
setMethod("drop_flagged", signature = c(object = "MetaboSet"),
  function(object, all_features = FALSE) {
    object <- .check_object(object, feature_flag = TRUE)
    if (!all_features) {
      object <- object[is.na(flag(object)), ]
    }
    if (!is.null(attr(object, "original_class"))) {
      object <- as(object, "MetaboSet")
      attr(object, "original_class") <- NULL
    }
    object
  }
)

#' @rdname drop_flagged
#' @export
setMethod("drop_flagged", signature = c(object = "SummarizedExperiment"),
  function(object, all_features = FALSE) {
    object <- .check_object(object)
    if (!all_features) {
      object <- object[is.na(flag(object)), ]
    }
    object
  }
)

#' Partially replace peak table with new values
#'
#' Replaces a subset of data in peak table of an object by new values.
#' Used after an operation such as imputation is computed only for a subset of 
#' features or samples such as only good-quality features that have not been 
#' flagged. This function is mainly used internally, but can be useful.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param y matrix containing new values to be merged into peak table
#'
#' @return A SummarizedExperiment or MetaboSet object with the new peak table 
#' values.
#'
#' @examples
#' data(example_set)
#' ex_set <- example_set[1:5, 1:5]
#' assay(ex_set)
#' # Create a matrix of replacment values for rows 1, 3, 5 and columns 1, 3, 4
#' replacement <- matrix(1:9,
#'   ncol = 3,
#'   dimnames = list(
#'     rownames(ex_set)[c(1, 3, 5)],
#'     colnames(ex_set)[c(1, 3, 4)]
#'   )
#' )
#' replacement
#' merged <- merge_assay(ex_set, replacement)
#' assay(merged)
#'
#' @noRd
merge_assay <- function(object, y, assay.type = NULL, name = NULL) {
  # Colnames and rownames should be found in the object
  if (!all(colnames(y) %in% colnames(assay(object))) || is.null(colnames(y))) {
    stop("Column names of y do not match column names of assay(object).")
  }
  if (!all(rownames(y) %in% rownames(assay(object))) || is.null(rownames(y))) {
    stop("Row names of y do not match row names of assay(object).")
  }
  from_to <- .get_from_to_names(object, assay.type, name)
  assay_tmp <- assay(object, from_to[[1]])
  assay_tmp[rownames(y), colnames(y)] <- y
  assay(object, from_to[[2]]) <- assay_tmp
  object
}

#' Impute missing values using random forest
#'
#' Impute the missing values in the peak table of the object using a
#' random forest. The estimated error in the imputation is logged.
#' It is recommended to set the seed number for reproducibility
#' (it is called random forest for a reason).
#' This a wrapper around \code{\link[missForest]{missForest}}.
#' Use parallelize = "variables" to run in parallel for faster testing.
#' NOTE: running in parallel prevents user from setting a seed number.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before imputation.
#' @param assay.type character, assay to be used in case of multiple assays
#' @param name character, name of the resultant assay in case of multiple assays
#' @param ... passed to \code{\link[missForest]{missForest}}
#'
#' @return An object as the one supplied, with missing 
#' values imputed.
#'
#' @examples
#' data(example_set)
#' missing <- mark_nas(example_set, 0)
#' set.seed(38)
#' imputed <- impute_rf(missing)
#'
#' @seealso \code{\link[missForest]{missForest}} for detail about the algorithm 
#' and the parameters
#'
#' @export
impute_rf <- function(object, all_features = FALSE, assay.type = NULL, 
                      name = NULL, ...) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("Package \'missForest\' needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("missForest package was used for random forest imputation:",
                citation("missForest"))
  # Start log
  log_text(paste("\nStarting random forest imputation at", Sys.time()))
  # Drop flagged features
  dropped <- drop_flagged(object, all_features)
  from_to <- .get_from_to_names(object, assay.type, name)
  dropped <- .check_object(dropped, assay.type = from_to[[1]])
  object <- .check_object(object, assay.type = from_to[[1]])

  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("missForest package not found.")
  }

  # Impute missing values
  mf <- missForest::missForest(xmis = t(assay(dropped, from_to[[1]])), ...)
  imputed <- t(mf$ximp)
  # Log imputation error
  log_text(paste0("Out-of-bag error in random forest imputation: ",
                  round(mf$OOBerror, digits = 3)))
  # Assign imputed data to the droppped
  rownames(imputed) <- rownames(assay(dropped))
  colnames(imputed) <- colnames(assay(dropped))
  # Attach imputed abundances to object
  object <- merge_assay(object, imputed, assay.type = from_to[[1]], 
                        name = from_to[[2]])
  log_text(paste("Random forest imputation finished at", Sys.time(), "\n"))
  
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}

#' Simple imputation
#'
#' Impute missing values using a simple imputation strategy. All missing values
#' of a feature are imputed with the same value. It is possible
#' to only impute features with a large number of missing values this way. This 
#' can be useful for using this function before random forest imputation to 
#' speed things up.
#' The imputation strategies available are:
#' \itemize{
#' \item a numeric value: impute all missing values in all features with the 
#' same value, e.g. 1
#' \item "mean": impute missing values of a feature with the mean of observed 
#' values of that feature
#' \item "median": impute missing values of a feature with the median of 
#' observed values of that feature
#' \item "min": impute missing values of a feature with the minimum observed 
#' value of that feature
#' \item "half_min": impute missing values of a feature with half the minimum 
#' observed value of that feature
#' \item "small_random": impute missing values of a feature with random numbers 
#' between 0 and the
#' minimum of that feature (uniform distribution, remember to set the seed 
#' number!).
#' }
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param value the value used for imputation, either a numeric or one of 
#' '"min", "half_min", "small_random", see above
#' @param na_limit only impute features with the proportion of NAs over this 
#' limit. For example, if \code{na_limit = 0.5}, only features with at least 
#' half of the values missing are imputed.
#' @param assay.type character, assay to be used in case of multiple assays
#' @param name character, name of the resultant assay in case of multiple assays
#'
#' @return A SummarizedExperiment or Metaboset object with imputed peak table.
#'
#' @examples
#' data(example_set)
#' missing <- mark_nas(example_set, 0)
#' imputed <- impute_simple(missing, value = "min")
#'
#' @export
impute_simple <- function(object, value, na_limit = 0, assay.type = NULL,
                          name = NULL) {
  from_to <- .get_from_to_names(object, assay.type, name)
  object <- .check_object(object, assay.type = from_to[[1]])
  imp <- assay(object, from_to[[1]])
  nas <- apply(imp, 1, prop_na)
  imp <- imp[nas > na_limit, , drop = FALSE]
  if (nrow(imp) == 0) {
    warning("None of the features satisfy the NA limit,", 
            " returning the original object.")
    return(object)
  }

  # Replace all missing values with the given constant
  if (is.numeric(value)) {
    imp[is.na(imp)] <- value
  } else if (value == "mean") {
    imp <- t(apply(imp, 1, function(x) {
      x[is.na(x)] <- finite_mean(x)
      x
    }))
  } else if (value == "median") {
    imp <- t(apply(imp, 1, function(x) {
      x[is.na(x)] <- finite_median(x)
      x
    }))
  } else if (value == "min") {
    imp <- t(apply(imp, 1, function(x) {
      x[is.na(x)] <- finite_min(x)
      x
    }))
  } else if (value == "half_min") {
    imp <- t(apply(imp, 1, function(x) {
      x[is.na(x)] <- finite_min(x) / 2
      x
    }))
  } else if (value == "small_random") {
    imp <- t(apply(imp, 1, function(x) {
      x[is.na(x)] <- stats::runif(n = sum(is.na(x)), 
                                  min = 0, max = finite_min(x))
      x
    }))
  } else {
    stop("value should be a numeric value or one of 'min',",
         " 'half_min', 'small_random'.")
  }

  obj <- merge_assay(object, imp, assay.type = from_to[[1]],
                     name = from_to[[2]])
  
  if (!is.null(attr(obj, "original_class"))) {
    object <- as(obj, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  obj
}

#' Inverse-rank normalization
#'
#' Applies inverse rank normalization to all features to approximate
#' a normal distribution.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param assay.type character, assay to be used in case of multiple assays
#' @param name character, name of the resultant assay in case of multiple assays
#' @return An object as the one supplied, with normalized features.
#'
#' @examples
#' data(example_set)
#' normalized <- inverse_normalize(example_set)
#' @export
inverse_normalize <- function(object, assay.type = NULL, name = NULL) {
  from_to <- .get_from_to_names(object, assay.type, name)
  object <- .check_object(object, assay.type = from_to[[1]])
  assay(object, from_to[[2]]) <- assay(object, from_to[[1]]) %>%
    apply(1, function(x) {
      stats::qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
    }) %>%
    t()
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}


#' A report of flagged features
#'
#' Computes the number of features at each stage of flagging for each mode.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#'
#' @return A data frame with the number of features at each stage of flagging.
#'
#' @examples
#' data(example_set)
#' flagged <- example_set %>%
#'   mark_nas(0) %>%
#'   flag_detection(group = "Group") %>%
#'   flag_quality()
#' flag_report(flagged)
#'
#' @export
flag_report <- function(object) {
  object <- .check_object(object, feature_split = TRUE)
  splits <- sort(unique(rowData(object)$Split))
  report <- data.frame()
  flag(object)[is.na(flag(object))] <- "Kept"
  for (split in splits) {
    tmp <- object[rowData(object)$Split == split, ]
    report_row <- flag(tmp) %>%
      table() %>%
      as.matrix() %>%
      t()
    report_row <- data.frame(Split = split, report_row)
    if (is.null(report_row$Kept)) {
      report_row$Kept <- 0
    }
    report_row$Total <- nrow(tmp)
    report_row$Flagged <- report_row$Total - report_row$Kept
    report <- dplyr::bind_rows(report, report_row)
  }
  report
}


# ---------- Logarithms ----------

#' Logarithm
#'
#' Log-transforms the peak table of SummarzedExperiment or MetaboSet object. 
#' Shortcuts for log2 and log10 also implemented.
#' For more information, see \code{\link{log}}.
#'
#' @param x a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param base the base of the logarithm
#'
#' @return A SummarizedExperiment or MetaboSet object with peak table
#' transformed.
#'
#' @name log
NULL

#' @rdname log
#' @export
setMethod("log", "MetaboSet", 
  function(x, base = exp(1)) {
    exprs(x) <- log(exprs(x), base = base)
    x
  }
)

#' @rdname log
#' @export
setMethod("log2", "MetaboSet", 
  function(x) {
    exprs(x) <- log2(exprs(x))
    x
  }
)

#' @rdname log
#' @export
setMethod("log10", "MetaboSet", 
  function(x) {
    exprs(x) <- log10(exprs(x))
    x
  }
)

#' Scale exprs data
#'
#' Applies the base R function scale to transposed peak table. 
#' See \code{\link[base]{scale}} for details.
#'
#' @param x a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param center,scale as in base scale function
#' 
#' @return A SummarizedExperiment or MetaboSet object with modified peak table.
#'
#' @name scale
#'
#' @export
NULL

#' @rdname scale
#' @export
setMethod("scale", "MetaboSet", 
  function(x, center = TRUE, scale = TRUE) {
    exprs(x) <- t(scale(t(exprs(x)), center = center, scale = scale))
    x
  }
)

#' Exponential function
#'
#' Apply exponential function to peak table.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param base base of the exponential
#'
#' @return An object with altered feature abundances.
#'
#' @examples 
#' data(example_set)
#' example_set <- mark_nas(example_set, value = 0)
#' log_data <- log(example_set)
#' orig_data <- exponential(log_data)
#'
#' @export
setGeneric("exponential", signature = "object",
           function(object, base = exp(1)) standardGeneric("exponential"))

#' @rdname exponential
#' @export
setMethod("exponential", c(object = "MetaboSet"),
  function(object, base = exp(1)) {
    exprs(object) <- base^exprs(object)
    object
  }
)


# ---------- Logarithms ----------

#' @rdname log
#' @export
setMethod("log", "SummarizedExperiment", 
  function(x, base = exp(1)) {
    assay(x) <- log(assay(x), base = base)
    x
  }
)

#' @rdname log
#' @export
setMethod("log2", "SummarizedExperiment", 
  function(x) {
    assay(x) <- log2(assay(x))
    x
  }
)

#' @rdname log
#' @export
setMethod("log10", "SummarizedExperiment", 
  function(x) {
    assay(x) <- log10(assay(x))
    x
  }
)

#' @rdname scale
#' @export
setMethod("scale", "SummarizedExperiment", 
  function(x, center = TRUE, scale = TRUE) {
    assay(x) <- t(scale(t(assay(x)), center = center, scale = scale))
    x
  }
)

#' @rdname exponential
#' @export
setMethod("exponential", c(object = "SummarizedExperiment"),
  function(object, base = exp(1)) {
    assay(object) <- base^assay(object)
    object
  }
)



#' Probabilistic quotient normalization
#'
#' Apply probabilistic quotient normalization (PQN) to the peak table of a 
#' SummarizedExperiment or MetaboSet object. By default, reference is 
#' calculated from high-quality QC samples and the median of the reference is 
#' used for normalization. Check parameters for more options.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param ref character, the type of reference samples to use for normalization.
#' @param method character, the method to use for calculating the reference 
#' sample.
#' @param all_features logical, should all features be used for calculating the 
#' reference sample?
#' @param assay.type character, assay to be used in case of multiple assays
#' @param name character, name of the resultant assay in case of multiple assays
#'
#' @return A SummarizedExperiment or MetaboSet object with altered feature 
#' abundances.
#'
#' @examples
#' data(example_set)
#' pqn_set <- pqn_normalization(example_set)
#'
#' @export
pqn_normalization <- function(object, ref = c("qc", "all"),
                              method = c("median", "mean"), 
                              all_features = FALSE, assay.type = NULL,
                              name = NULL) {
  log_text("Starting PQN normalization")
  ref <- match.arg(ref)
  method <- match.arg(method)
  # Use only good-quality features for calculating reference spectra
  from_to <- .get_from_to_names(object, assay.type, name)
  object <- .check_object(object, pheno_QC = TRUE, assay.type = from_to[[1]])
  ref_data <- assay(drop_flagged(object, all_features), from_to[[1]])
  
  
  # Select reference samples
  switch(ref, qc = reference <- ref_data[, object$QC == "QC"],
         all = reference <- ref_data)
  if (ncol(reference) == 0 || all(is.na(reference))) {
    stop("No specified reference samples found.")
  }
  # Calculate reference spectrum
  switch(method,
         median = reference_spectrum <- apply(reference, 1, finite_median),
         mean = reference_spectrum <- apply(reference, 1, finite_mean))
  log_text(paste("Using", method, "of", ref, "samples as reference spectrum"))
  # Calculate median of quotients
  quotients <- ref_data / reference_spectrum
  quotient_md <- apply(quotients, 2, finite_median)
  # Do the normalization
  data <- assay(object, from_to[[1]])
  pqn_data <- t(t(data) / quotient_md)
  colnames(pqn_data) <- colnames(data)
  rownames(pqn_data) <- rownames(data)
  assay(object, from_to[[2]]) <- pqn_data
  
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}
