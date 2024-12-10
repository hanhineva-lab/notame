#' Save plot to chosen format
#'
#' Saves the given plot to a file. Supports pdf, svg, emf, png and tiff formats.
#' If an error occurs with the plot, an empty file is created.
#'
#' @param p a ggplot object
#' @param file the file path
#' @param ... other arguments to plot function, like width and height
#'
#' @return None, the function is invoked for its plot-saving side effect.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' lm_results <- perform_lm(drop_qcs(example_set), 
#'   formula_char = "Feature ~ Group")
#'
#' p <- volcano_plot(lm_results,
#'   x = "GroupB_Estimate",
#'   p = "GroupB_P", p_fdr = "GroupB_P_FDR",
#'   label = "Feature_ID",
#'   fdr_limit = 0.1
#' )
#'
#' save_plot(p, file = "test.pdf")
#' \dontshow{setwd(.old_wd)}
#'
#' @seealso \code{\link[grDevices]{pdf}},
#' \code{\link[devEMF]{emf}},
#' \code{\link[grDevices]{svg}},
#' \code{\link[grDevices]{png}},
#' \code{\link[grDevices]{tiff}}
#'
#' @export
save_plot <- function(p, file, ...) {
  # Create folder automatically
  folder <- dirname(file)
  if (!file.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }

  format <- utils::tail(unlist(strsplit(basename(file), split = "\\.")), 1)
  switch(format,
    "emf" = devEMF::emf(file, ...),
    "pdf" = grDevices::pdf(file, ...),
    "svg" = grDevices::svg(file, ...),
    "png" = grDevices::png(file, ...),
    "tiff" = grDevices::tiff(file, ...),
    stop("File format '", format, "' is not valid, saving failed"))
  tryCatch(
    {
      plot(p)
      grDevices::dev.off()
      log_text(paste("Saved to:", file))
    },
    error = function(e) {
      grDevices::dev.off()
      stop(e$message, call. = FALSE)
    }
  )
}

# Helper function for handling errors and keeping track of file names
.save_name <- function(object, prefix, format, fun, name, 
                       width = 7, height = 7, ...) {
  save_seed <- .Random.seed
  file_names <- ""
  p <- NULL
  tryCatch(
    {
      p <- fun(object, ...)
    },
    error = function(e) {
      message("Problem with plot named ", name, ":\n", e$message)
    }
  )
  
  if (!is.null(p)) {
    file_name <- paste0(prefix, "_", name, ".", format)
    save_plot(p, file = file_name, width = width, height = height)
    assign("file_names", paste(file_names, file_name))
  }
  assign(".Random.seed", save_seed, envir = .GlobalEnv)

}

.merge_to_pdf <- function(prefix, file_names, remove_singles) {
  prefix <- gsub("_$", "", prefix)
  merged_file <- paste0(prefix, ".pdf")
  os <- Sys.info()[["sysname"]]
  output <- NULL
  if (os == "Windows") {
    # Merge files
    output <- shell(paste("pdftk", file_names, "cat output", merged_file),
                    intern = TRUE)
  } else if (os == "Linux") {
    output <- system(paste("pdfunite", file_names, merged_file),
                     intern = TRUE)
  } else if (os == "Darwin") {
    output <- system(paste('"/System/Library/Automator/Combine PDF',
                           'Pages.action/Contents/Resources/join.py" -o',
                           merged_file, file_names),
                     intern = TRUE)
  } else {
    log_text(paste0("Unfortunately your operating system is",
                    " not yet supported by the merging"))
    return()
  }
  if (length(output) && output != "0") {
    log_text(paste("Merging plots resulted in the following message:",
                   paste0(output, collapse = " ")))
  } else {
    log_text(paste("Attempted merging plots to", merged_file))
    if (remove_singles) {
      log_text("Removing single plot files")
      if (os == "Windows") {
        output2 <- shell(paste("del", file_names), intern = TRUE)
      } else {
        output2 <- system(paste("rm", file_names), intern = TRUE)
      }
      if (length(output2) && output2 != "0") {
        log_text(paste("Removing single plot files resulted in",
                       "the following message:",
                        paste0(output2, collapse = " ")))
      }
    }
  }
}


#' Write all relevant pretreatment visualizations to pdf
#'
#' A wrapper around all the major visualization functions, used for visualizing 
#' data between major steps of data preprocessing. Saves all visualizations as 
#' PDFs with a set prefix on filenames.
#'
#' @param object A SummarizedExperiment or MetaboSet object
#' @param prefix character, a file path prefix added to the file paths
#' @param format character, format in which the plots should be saved, DOES NOT 
#' support raster formats
#' @param perplexity perplexity for t-SNE plots
#' @param merge logical, whether the files should be merged to a single PDF, 
#' see Details
#' @param remove_singles logical, whether to remove single plot files after 
#' merging. Only used if \code{merge = TRUE}
#' @param group character, name of pheno data column containing the group labels
#' @param time character, name of pheno data column containing timepoints
#' @param id character, name of pheno data column containing subject identifiers
#' @param color character, name of pheno data column used for coloring sample
#' labels for dendrograms
#' @param assay.type character, assay to be used in case of multiple assays
#' @return None, the function is invoked for its plot-saving side effect.
#'
#' @details If \code{merge} is \code{TRUE} and \code{format} is \code{pdf},
#' then a file containing all the visualizations named \code{prefix.pdf} will 
#' be created. NOTE: on Windows this requires installation of pdftk 
#' (\url{https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/})
#' and on Linux you need to have pdfunite installed. On MacOS, no external 
#' software is needed. Note that at least on Windows, prefix should be a path 
#' from the root, so that the underlying system command will find the files.
#' The type of visualizations to be saved depends on the type of object.
#' Here is a comprehensive list of the visualizations:
#' \itemize{
#' \item Distribution of quality metrics and flags \code{\link{plot_quality}}
#' \item Boxplots of each sample in injection order 
#' \code{\link{plot_sample_boxplots}}
#' \item PCA scores plot of samples colored by injection order 
#' \code{\link{plot_pca}}
#' \item t-SNE plot of samples colored by injection order 
#' \code{\link{plot_tsne}}
#' \item If the object has over 60 samples, hexbin versions of the PCA and t-
#' SNE plots above
#' \code{\link{plot_pca_hexbin}}, \code{\link{plot_tsne_hexbin}}
#' \item Dendrogram of samples ordered by hierarchical clustering, sample 
#' labels colored by group if present
#' \code{\link{plot_dendrogram}}
#' \item heat map of intersample distances, ordered by hierarchical clustering 
#' \code{\link{plot_sample_heatmap}}
#' \item If the object has QC samples: \itemize{
#' \item Density function of the intersample distances in both QCs and 
#' biological samples \code{\link{plot_dist_density}}
#' \item Histograms of p-values from linear regression of features against 
#' injection order in both QCs and biological samples 
#' \code{\link{plot_p_histogram}}}
#' \item If the object has a group column: \itemize{
#' \item PCA and tSNE plots with points shaped and colored by group 
#' \code{\link{plot_pca}}, \code{\link{plot_tsne}}
#' }
#' \item If the object has a time column: \itemize{
#' \item PCA and tSNE plots with points shaped and colored by time 
#' '\code{\link{plot_pca}}, \code{\link{plot_tsne}}
#' \item Dendrogram of samples ordered by hierarchical clustering, sample 
#' labels colored by time point \code{\link{plot_dendrogram}}
#' }
#' \item If the object has a group column OR a time column: \itemize{
#' \item Boxplots of samples ordered and colored by group and/or time 
#' \code{\link{plot_sample_boxplots}}
#' }
#' \item If the object has a group column AND a time column: \itemize{
#' \item PCA and tSNE plots with points shaped by group and colored by time
#' \code{\link{plot_pca}}, \code{\link{plot_tsne}}
#' }
#' \item If the object has a time column AND a subject column: \itemize{
#' \item PCA and tSNE plots with arrows connecting the samples of each subject 
#' in time point order
#' \code{\link{plot_pca_arrows}}, \code{\link{plot_tsne_arrows}}
#' }
#' }
#'
#' @seealso \code{\link[notame]{save_plot}}
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' visualizations(example_set, prefix="figures/example_set", perplexity=5,
#'                group = "Group", color = "Group", time = "Time", 
#'                id = "Subject_ID")
#' \dontshow{setwd(.old_wd)}
#'
#' @export
visualizations <- function(object, prefix, format = "pdf", perplexity = 30,
                           merge = FALSE, remove_singles = FALSE, group = NULL, 
                           time = NULL, id = NULL, color = NULL,
                           assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_QC = TRUE,
                         pheno_cols = c(time, id, color), assay.type = from)
  assays(object) <- assays(object)[from]
  file_names <- ""
  if (sum(object$QC == "QC")) {
    .save_name(object, prefix, format, fun = plot_dist_density, 
               name = "density_plot", width = 8, height = 6)
    .save_name(object, prefix, format, plot_injection_lm, "lm_p_histograms")
  }
  # Quality metrics
  .save_name(object, prefix, format, plot_quality, "quality_metrics")
  # Plots with injection order
  .save_name(object, prefix, format, plot_sample_boxplots,
             "boxplots_injection", order_by = "Injection_order", 
             fill_by = "QC", width = 15)
  .save_name(object, prefix, format, plot_pca, "PCA_injection", 
             color = "Injection_order")
  .save_name(object, prefix, format, plot_tsne, "tSNE_injection",
             perplexity = perplexity, color = "Injection_order")
  # Clustering
  .save_name(object, prefix, format, plot_dendrogram, "dendrogram", 
             width = 15, color = color)
  .save_name(object, prefix, format, plot_sample_heatmap, 
             "heatmap_samples", width = 15, height = 16, group = group)
  # For large sets, plot hexbin plots
  if (ncol(object) > 60) {
    .save_name(object, prefix, format, plot_pca_hexbin, "PCA_hexbin")
    .save_name(object, prefix, format, plot_tsne_hexbin, 
               "tSNE_hexbin", perplexity = perplexity)
  }
  # If not grouped, plot PCA and t-SNE on QC information
  if (is.null(colData(object)[, group])) {
    group <- "QC"
  }
  .save_name(object, prefix, format, plot_pca, "PCA_group", color = group)
  .save_name(object, prefix, format, plot_tsne, 
             "tSNE_group", perplexity = perplexity, color = group)
  # Time point
  if (!is.null(colData(object)[, time])) {
    .save_name(object, prefix, format, plot_pca, "PCA_time",
               color = time)
    .save_name(object, prefix, format, plot_tsne, "tSNE_time",
               color = time, perplexity = perplexity)
    .save_name(object, prefix, format, plot_dendrogram, "dendrogram_time",
               color = time, width = 15)
  }
  # Time point OR group
  if (!is.null(colData(object)[, group]) || !is.null(colData(object)[, time])){
    by <- c(group, time)
    .save_name(object, prefix, format, plot_sample_boxplots, 
               "boxplots_group", width = 15, 
               order_by = by, fill_by = by)
  }
  # Time point AND group
  if (!is.null(colData(object)[, group]) && !is.null(colData(object)[, time])) {
    .save_name(object, prefix, format, plot_pca, "PCA_group_time",
               color = time, shape = group)
    .save_name(object, prefix, format, plot_tsne, "tSNE_group_time", 
               color = time, shape = group,
               perplexity = perplexity)
  }
  # Multiple time points per subject
  if (!is.null(colData(object)[, time]) &&
    !is.null(colData(object)[, id]) &&
    sum(object$QC == "QC") == 0) {
    .save_name(object, prefix, format, plot_pca_arrows, "PCA_arrows", 
               color = group, time = time, subject = id)
    .save_name(object, prefix, format, plot_tsne_arrows, "tSNE_arrows", 
               perplexity = perplexity, color = group, time = time, 
               subject = id)
  }
  if (merge && format == "pdf") {
    .merge_to_pdf(prefix, file_names, remove_singles)
  }
}
