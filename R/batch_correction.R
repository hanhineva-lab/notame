#' Remove Unwanted Variation (RUV) between batches
#'
#' An interface for \code{link[RUVSeq]{RUVs}} method
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param replicates list of numeric vectors, indexes of replicates
#' @param k The number of factors of unwanted variation to be estimated from 
#' the data.
#' @param assay.type character, assay to be used in case of multiple assays
#' @param name character, name of the resultant assay in case of multiple assays
#' @param ... other parameters passed to \code{link[RUVSeq]{RUVs}}
#' @return A SummarizedExperiment or Metaboset object with the normalized data.
#'
#' @examples
#' data(example_set)
#' # Batch correction
#' replicates <- list(which(example_set$QC == "QC"))
#' batch_corrected <- ruvs_qc(example_set, replicates = replicates)
#' # Evaluate batch correction
#' pca_bhattacharyya_dist(example_set, batch = "Batch")
#' pca_bhattacharyya_dist(batch_corrected, batch = "Batch")
#'
#' @export
ruvs_qc <- function(object, replicates, k = 3, 
                    assay.type = NULL, name = NULL, ...) {
  if (!requireNamespace("RUVSeq", quietly = TRUE)) {
    stop("Bioconductor package RUVSeq needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  
  .add_citation("RUVSeq was used for batch correction:", citation("RUVSeq"))    
  
  from_to <- .get_from_to_names(object, assay.type, name)
  object <- .check_object(object, assay.type = from_to[[1]])
  
  # Transform data to pseudo counts for RUVs
  assay(object, from_to[[1]])[assay(object, from_to[[1]]) == 0] <- 1
  assay(object) <- round(assay(object))
  # Pad each replicate vector with -1 and transform to matrix
  max_len <- max(vapply(replicates, length, integer(1)))
  scIdx <- matrix(-1, nrow = length(replicates), ncol = max_len)
  # Populate matrix of replicates
  for (i in seq_along(replicates)) {
    scIdx[i, seq_along(replicates[[i]])] <- replicates[[i]]
  }
  # Perform batch correction
  ruv_results <- RUVSeq::RUVs(x = assay(object), cIdx = rownames(object),
                              k = k, scIdx = scIdx, ...)
  # Include results in object
  assay(object, from_to[[2]]) <- ruv_results$normalizedCounts
  colData(object) <- cbind(colData(object), ruv_results$W)
  
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}


#' Bhattacharyya distance between batches in PCA space
#'
#' Computes Bhattacharyya distance between all pairs of batches using 
#' \code{\link[fpc:bhattacharyya.matrix]{fpc:bhattacharyya.matrix}} 
#' after projecting the samples into PCA space with 
#' \code{\link[pcaMethods]{pca}}.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param batch column name of pheno data giving the batch labels
#' @param all_features logical, should all features be used? If FALSE
#' (the default), flagged features are removed before imputation.
#' @param center logical, should the data be centered prior to PCA?
#' (usually yes)
#' @param scale scaling used, as in 
#' \code{\link[pcaMethods]{prep}}. Default is "uv" for unit variance
#' @param nPcs the number of principal components to use
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters to \code{\link[pcaMethods]{pca}}
#'
#' @return A matrix of Bhattacharyya distances between batches.
#'
#' @examples
#' data(example_set)
#' # Batch correction
#' replicates <- list(which(example_set$QC == "QC"))
#' batch_corrected <- ruvs_qc(example_set, replicates = replicates)
#' # Evaluate batch correction
#' pca_bhattacharyya_dist(example_set, batch = "Batch")
#' pca_bhattacharyya_dist(batch_corrected, batch = "Batch")
#'
#' @export
pca_bhattacharyya_dist <- function(object, batch, all_features = FALSE, 
                                   center = TRUE, scale = "uv", nPcs = 3, 
                                   assay.type = NULL, ...){ 
  if (!requireNamespace("fpc", quietly = TRUE)) {
    stop("Package \"fpc\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \"pcaMethods\" needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  .add_citation("PCA was performed using pcaMethods package:",
                citation("pcaMethods"))
  .add_citation("fpc package was used for Bhattacharyaa distance computation:",
                citation("fpc"))
  from <- .get_from_name(object, assay.type)
  # Drop flagged features if not told otherwise
  object <- .check_object(object, pheno_factors = batch, assay.type = from)
  object <- drop_flagged(object, all_features)
  # PCA to 2 dimensions
  pca_res <- pcaMethods::pca(t(assay(object, from)), center = center, 
                             scale = scale, nPcs = nPcs, ...)
  pca_scores <- pcaMethods::scores(pca_res)
  # Split to batches
  batches <- list()
  for (b in unique(colData(object)[, batch])) {
    batches[[b]] <- pca_scores[colData(object)[, batch] == b, ]
  }
  # Compute means and covariance matrices for Bhattacharyya distance
  muarray <- vapply(batches, colMeans, double(nPcs))
  sigmaarray <- array(vapply(batches, stats::cov, double(nPcs * nPcs)), 
                      dim = c(nPcs, nPcs, length(batches)))
  # Compute Bhattacharyya distance
  fpc::bhattacharyya.matrix(muarray, sigmaarray, ipairs = "all",
                            misclassification.bound = FALSE)
}


.pooled_variance <- function(x, group) {
  # Remove missing values
  group <- group[!is.na(x)]
  x <- x[!is.na(x)]
  # Split to groups
  group_list <- split(x, group)
  n_1 <- vapply(group_list, length, integer(1)) - 1 # n - 1
  # Pooled variance
  sum(n_1 * vapply(group_list, var, numeric(1))) / sum(n_1)
}

.between_variance <- function(x, group) {
  # Remove missing values
  group <- group[!is.na(x)]
  x <- x[!is.na(x)]
  # Split to groups
  group_list <- split(x, group)
  n <- vapply(group_list, length, integer(1))
  means <- vapply(group_list, mean, numeric(1))
  k_1 <- length(unique(group)) - 1
  # Between group variance formula
  sum(n * (means - mean(x))^2) / k_1
}

.repeatability <- function(x, group) {
  # Calculate pooled variance
  pv <- .pooled_variance(x, group)
  # Calculate between group variance
  bv <- .between_variance(x, group)
  bv / (bv + pv)
}


#' Compute repeatability measures
#'
#' Computes repeatability for each feature with the following formula:
#' \deqn{\frac{\sigma^2_{between}}{\sigma^2_{between} + \sigma^2_{within}}}
#' The repeatability ranges from 0 to 1. Higher repeatability depicts less
#' variation between batches.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param group column name of pheno data giving the group labels
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A data frame with one row per feature with the repeatability measure.
#'
#' @examples
#' data(example_set)
#' # Batch correction
#' replicates <- list(which(example_set$QC == "QC"))
#' batch_corrected <- ruvs_qc(example_set, replicates = replicates)
#' # Evaluate batch correction
#' rep_orig <- perform_repeatability(example_set, group = "Group")
#' mean(rep_orig$Repeatability, na.rm = TRUE)
#' rep_corr <- perform_repeatability(batch_corrected, group = "Group")
#' mean(rep_corr$Repeatability, na.rm = TRUE)
#'
#' @export
perform_repeatability <- function(object, group, assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_factors = group, assay.type = from)
  group <- colData(object)[, group]
  features <- rownames(object)
  repeatability <- BiocParallel::bplapply(
    X = features, 
    FUN = function(feature) {
      result_row <- data.frame(
        Feature_ID = feature,
        Repeatability = .repeatability(assay(object)[feature, ], group))
      }
    )
  do.call(rbind, repeatability)
}

#' Save batch correction plots
#'
#' Saves plots of each feature showing the effect of batch correction.
#' Plots show QC samples and regular samples inside each batch, plus the
#' batch mean for biological samples and QC samples as a horizontal line.
#' The dashed line represents QC mean, the filled line represents biological
#' sample mean.
#' NOTE: if you change the shape variable, be sure to set a shape scale as well,
#' the default scale only has 2 values, so it can only accomodate 2 shapes.
#'
#' @param orig,corrected \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} objects before and after batch effect correction
#' @param file path to the PDF file where the plots will be saved
#' @param width,height width and height of the plots in inches
#' @param batch,color,shape column names of pheno data for batch labels,
#' and column used for coloring and shaping points (by default batch and QC)
#' @param color_scale,shape_scale scales for color and scale as returned by 
#' ggplot functions.
#' @param assay.type1 character, assay of orig to be used in case of 
#' multiple assays.
#' @param assay.type2 character, assay of corrected to be used in case of
#' multiple assays. If corrected is not supplied, this argument selects
#' another assay from orig.
#' @return None, the function is invoked for its plot-saving side effect.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' data(example_set)
#' # Batch correction
#' replicates <- list(which(example_set$QC == "QC"))
#' batch_corrected <- ruvs_qc(example_set, replicates = replicates)
#' # Plots of each feature
#' save_batch_plots(
#'   orig = example_set[1:10], corrected = batch_corrected[1:10],
#'   file = "batch_plots.pdf"
#' )
#' \dontshow{setwd(.old_wd)}
#' @export
save_batch_plots <- function(orig, corrected, file, width = 14, 
                             height = 10, batch = "Batch", color = "Batch", 
                             shape = "QC",
                             color_scale = getOption("notame.color_scale_dis"),
                             shape_scale = 
                             scale_shape_manual(values = c(15, 21)),
                             assay.type1 = NULL, assay.type2 = NULL) {
                               
  if (missing(corrected) && is(orig, "SummarizedExperiment")) {
    from1 <- .get_from_name(orig, assay.type1)
    from2 <- .get_from_name(orig, assay.type2)
    orig <- .check_object(orig, pheno_factors = c(batch, shape),
                          pheno_cols = color)
    data_orig <- combined_data(orig, from1)
    data_corr <-  combined_data(orig, from2)
  } else {
    from1 <- .get_from_name(orig, assay.type1)
    from2 <- .get_from_name(corrected, assay.type2)
    orig <- .check_object(orig, pheno_factors = c(batch, shape),
                          pheno_cols = color)
    corrected <- .check_object(corrected, pheno_factors = c(batch, shape),
                               pheno_cols = color)
    data_orig <- combined_data(orig, assay.type = from1)
    data_corr <- combined_data(corrected, from2)
  }
  
  # Prepare data.frame for batch means with batch and injection order range
  batch_injections <- data_orig %>%
    dplyr::group_by(!!dplyr::sym(batch)) %>%
    dplyr::summarise(start = min(.data$Injection_order), 
                     end = max(.data$Injection_order))

  batch_mean_helper <- function(data) {
    data %>%
      dplyr::group_by(!!dplyr::sym(batch)) %>%
      dplyr::summarise_at(rownames(orig), finite_mean) %>%
      dplyr::left_join(batch_injections, ., by = batch)
  }
  
  get_batch_means <- function(data) {
    batch_means <- batch_mean_helper(data) %>%
      dplyr::mutate(QC = "Sample")
    batch_means_qc <- data %>%
      dplyr::filter(.data$QC == "QC") %>%
      batch_mean_helper() %>%
      dplyr::mutate(QC = "QC")
    rbind(batch_means, batch_means_qc)
  }
  # Get batch means for QC and biological samples of the original data
  batch_means_orig <- get_batch_means(data_orig)
  # Get batch means for QC and biological samples of the corrected data
  batch_means_corr <- get_batch_means(data_corr)
  
  batch_plot_helper <- function(data, fname, batch_means) {
    p <- ggplot() +
      geom_point(data = data, 
                 mapping = aes(x = .data[["Injection_order"]], 
                               y = .data[[fname]],
                               color = .data[[color]], 
                               shape = .data[[shape]])) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      color_scale +
      shape_scale
    p <- p +
      geom_segment(data = batch_means, 
                   mapping = aes(x = .data[["start"]], xend = .data[["end"]],
                                 y = .data[[fname]], yend = .data[[fname]],
                                 color = .data[[color]], 
                                 linetype = .data[["QC"]]),
                             size = 1) +
      scale_linetype(guide = "none")
      
    p
  }
  # Save plots with original and corrected data to pdf
  grDevices::pdf(file, width = width, height = height)
  for (feature in rownames(orig)) {
    p1 <- batch_plot_helper(data_orig, feature, batch_means_orig)
    p2 <- batch_plot_helper(data_corr, feature, batch_means_corr)
    p <- cowplot::plot_grid(p1, p2, nrow = 2)
    plot(p)
  }
  grDevices::dev.off()
}
