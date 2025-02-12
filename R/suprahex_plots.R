#' Sample suprahex plots
#'
#' Plots supraHex plots of each sample using functions from the supraHex 
#' package. See the supraHex paper and package vignette for more information.
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param all_features if FALSE, flagged features are droppped
#' @param sample_labels the column for labels of samples in the plot
#' @param grid_xdim,grid_ydim dimensions of the grid for the samples
#' @param title.xy position of the sample label relative to the supraHex
#' @param title.rotate rotation of the sample label in degrees
#' @param height a numeric value specifying the height of device
#' @param fontsize the fontsize for sample labels
#' @param colormap colormap for the hexagons
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters for supraHex::sPipeline
#'
#' @return None, the function is invoked for its plotting side effect.
#'
#' @examples
#' data(example_set)
#' plot_sample_suprahex(example_set[, 1:20],
#'   xdim = 5, title.xy = c(0.35, 1),
#'   grid_xdim = 7, grid_ydim = 7, sample_labels = "Group"
#' )
#'
#' @seealso \code{\link[supraHex]{sPipeline}},
#' \code{\link[supraHex]{sCompReorder}},
#'  \code{\link[supraHex]{visCompReorder}}
#'
#' @export
plot_sample_suprahex <- function(object, all_features = FALSE,
                                 sample_labels = "Sample_ID", grid_xdim = NULL,
                                 grid_ydim = NULL, title.xy = c(0.35, 1), 
                                 title.rotate = 0, height = 7, fontsize = 10,
                                 colormap = "jet", assay.type = NULL, ...) {
  if (!requireNamespace("supraHex", quietly = TRUE)) {
    stop("Package \'supraHex\' needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("supraHex package was used for suprahexagonal maps:",
                citation("supraHex"))
  from <- .get_from_name(object, assay.type)
  object <- drop_flagged(object, all_features = all_features)
  object <- .check_object(object, pheno_cols = sample_labels, assay.type = from)
  
  data <- scale(assay(object, from))
  colnames(data) <- colData(object)[, sample_labels]

  s_map <- supraHex::sPipeline(data = data, ...)
  s_reorder <- supraHex::sCompReorder(sMap = s_map, xdim = grid_xdim, 
                                      ydim = grid_ydim)

  supraHex::visCompReorder(sMap = s_map, sReorder = s_reorder, newpage = FALSE,
                           colormap = colormap, title.xy = title.xy, 
                           title.rotate = title.rotate, height = height,
                           gp = grid::gpar(fontsize = fontsize))
}
