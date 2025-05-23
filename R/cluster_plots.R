#' Sample dendrogram
#'
#' Draws a dendrogram of a hierarchical clustering applied to the samples of an 
#' experiment.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param color character, name of the column used for coloring the sample 
#' labels
#' @param dist_method distance method used in clustering as in
#' \code{\link[stats]{dist}}
#' @param clust_method method used in clustering as in  
#' \code{\link[stats]{hclust}}
#' @param center logical, should the data be centered?
#' @param scale scaling used, as in 
#' \code{\link[pcaMethods]{prep}}. Default is "uv" for unit variance
#' @param title The plot title
#' @param subtitle The plot subtitle
#' @param color_scale the color scale as returned by a ggplot function.
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A ggplot object.
#'
#' @examples
#' data(example_set)
#' plot_dendrogram(example_set, color = "Group")
#'
#' @seealso \code{\link{dist}} \code{\link{hclust}}
#'
#' @export
plot_dendrogram <- function(object, all_features = FALSE, 
                            color, dist_method = "euclidean",
                            clust_method = "ward.D2",
                            center = TRUE, scale = "uv", 
                            title = "Dendrogram of hierarchical clustering",
                            subtitle = NULL, 
                            color_scale = getOption("notame.color_scale_dis"),
                            assay.type = NULL) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \'pcaMethods\' needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = color, assay.type = from)
  
  subtitle <- subtitle %||% paste("Distance method:", dist_method, 
                                  "Clustering method:", clust_method)
  # change name to matrix
  assay <- pcaMethods::prep(t(assay(object, from)), center = center, 
                            scale = scale)

  d_data <- stats::dist(assay, method = dist_method) %>%
    stats::hclust(method = clust_method) %>%
    stats::as.dendrogram() %>%
    ggdendro::dendro_data()

  labels <- ggdendro::label(d_data) %>%
    dplyr::mutate(label = .data$label) %>%
    dplyr::left_join(colData(object)[c("Sample_ID", color)], 
                     by = c("label" = "Sample_ID"), copy = TRUE)
  labels[, color] <- as.factor(labels[, color])
  p <- ggplot(ggdendro::segment(d_data)) +
    geom_segment(aes(x = .data$x, y = .data$y,
                     xend = .data$xend, yend = .data$yend)) +
    geom_text(data = labels, aes(x = .data[["x"]], y = .data[["y"]],
                                 label = .data[["label"]],
                                 color = .data[[color]]),
              angle = 90, hjust = 1) +
    ggdendro::theme_dendro() +
    color_scale +
    labs(title = title, subtitle = subtitle)

  p
}
#' Sample heatmap
#'
#' Draws a heatmap of the distances between the samples of an experiment,
#' the samples are ordered by hierarchical clustering.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param dist_method distance method used in clustering as in
#' \code{\link[stats]{dist}}
#' @param clust_method method used in clustering as in  
#' \code{\link[stats]{hclust}}
#' @param center logical, should the data be centered?
#' @param scale scaling used, as in 
#' \code{\link[pcaMethods]{prep}}. Default is "uv" for unit variance
#' @param group_bar logical, should a bar showing the groups be drawn under the 
#' heat map?
#' @param group character, name of the column used for coloring the group bar
#' @param title The plot title
#' @param subtitle The plot subtitle
#' @param fill_scale_con Continuous fill scale for the heatmap as returned by a 
#' ggplot function
#' @param fill_scale_dis Discrete fill scale for the group bar as returned by a 
#' ggplot function
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A ggplot object. If \code{group_bar} is \code{TRUE}, the plot will 
#' consist of multiple parts and is harder to modify.
#'
#' @examples
#' data(example_set)
#' plot_sample_heatmap(example_set, group = "Group")
#'
#' @seealso \code{\link{dist}} \code{\link{hclust}}
#'
#' @export
plot_sample_heatmap <- function(object, all_features = FALSE, 
                                dist_method = "euclidean", 
                                clust_method = "ward.D2",
                                center = TRUE, scale = "uv",
                                group_bar = TRUE, group = NULL,
                                title = "Heatmap of distances between samples",
                                subtitle = NULL, fill_scale_con =
                                getOption("notame.fill_scale_con"),
                                fill_scale_dis =
                                getOption("notame.fill_scale_dis"),
                                assay.type = NULL) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \'pcaMethods\' needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = group, assay.type = from)

  # Default settings
  subtitle <- subtitle %||% paste("Distance method:", dist_method, 
                                  "Clustering method:", clust_method)

  assay <- pcaMethods::prep(t(assay(object, from)), center = center,
                            scale = scale)

  # Distances
  distances <- stats::dist(assay, method = dist_method)
  # Hierarchical clustering for ordering
  hc <- stats::hclust(distances, method = clust_method)
  hc_order <- hc$labels[hc$order]

  # From wide to long format for ggplot
  distances_df <- as.matrix(distances) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("X") %>%
    tidyr::gather("Y", "Distance", -"X")
  # Set the correct order given by hclust
  distances_df$X <- factor(distances_df$X, levels = hc_order, ordered = TRUE)
  distances_df$Y <- factor(distances_df$Y, levels = rev(hc_order), 
                           ordered = TRUE)

  # Heatmap
  p <- ggplot(distances_df, aes(.data$X, .data$Y, fill = .data$Distance)) +
    geom_tile(color = NA) +
    fill_scale_con +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.05)) +
    coord_fixed()
  # Group bar
  if (group_bar && !is.null(group)) {
    pheno_data <- colData(object)
    pheno_data$Sample_ID <- factor(pheno_data$Sample_ID, levels = hc_order,
                                   ordered = TRUE)

    gb <- ggplot(pheno_data, aes(x = .data[["Sample_ID"]], y = 1, 
                                 fill = .data[[group]])) +
      geom_tile(color = "white") +
      theme_void() +
      fill_scale_dis

    p <- cowplot::plot_grid(p, gb, ncol = 1, align = "v", 
                            rel_heights = c(10 / 11, 1 / 11))
  }

  p
}
