# ------- HELPER FUNCTIONS ----------------

# Helper function for computing PCA
.pca_helper <- function(object, pcs, center, scale, assay.type, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \'pcaMethods\' needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("PCA was performed using pcaMethods package:",
                citation("pcaMethods"))
  res_pca <- pcaMethods::pca(t(assay(object, assay.type)),
                             nPcs = max(pcs), scale = scale, 
                             center = center, ...)
  pca_scores <- as.data.frame(pcaMethods::scores(res_pca))[, pcs]
  r2 <- summary(res_pca)["R2", pcs]
  labels <- paste0(paste0("PC", pcs), " (", scales::percent(r2), ")")

  return(list(pca_scores = pca_scores, labels = labels))
}

# Helper function for computing t-SNE
.t_sne_helper <- function(object, center, scale, perplexity, 
                          pca_method, assay.type, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \'pcaMethods\' needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Package \'Rtsne\' needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  .add_citation("Rtsne package was used for t-SNE figures:", citation("Rtsne"))

  prepd <- pcaMethods::prep(t(assay(object, assay.type)), 
                            center = center, scale = scale)

  if (sum(is.na(prepd)) > 0) {
    res_pca <- pcaMethods::pca(t(assay(object, assay.type)), 
                               method = pca_method, 
                               nPcs = min(nrow(object), ncol(object), 50),
                               scale = "none", center = FALSE)
    pca_scores <- pcaMethods::scores(res_pca)
    res_tsne <- Rtsne::Rtsne(pca_scores, perplexity = perplexity,
                             pca = FALSE, ...)
  } else {
    res_tsne <- Rtsne::Rtsne(prepd, perplexity = perplexity, ...)
  }
  data.frame(res_tsne$Y)
}

# -------------- SCATTER PLOTS ---------------

#' PCA scatter plot
#'
#' Computes PCA using one of the methods provided in the Bioconductor package
#' pcaMethods and plots the two first principal components.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in \code{\link[pcaMethods]{prep}}. Default is 
#' "uv" for unit variance
#' @param color character, name of the column used for coloring the points. Set 
#' to NULL for black color.
#' @param shape character, name of the column used for shape. Set to NULL for 
#' uniform round shapes.
#' @param point_size numeric, size of the points.
#' @param label character, name of the column used for point labels
#' @param density logical, whether to include density plots to both axes.
#' The density curves will be split and colored by the 'color' variable.
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function. Set to 
#' NA to choose the appropriate scale based on the class of the coloring 
#' variable.
#' @param shape_scale the shape scale as returned by a ggplot function
#' @param fill_scale the fill scale used for density curves.
#' If a continuous variable is used as color, density curve will be colorless.
#' @param text_base_size numeric, base size for text
#' @param point_size numeric, size of the points
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional arguments passed to \code{\link[pcaMethods]{pca}}
#'
#' @return A ggplot object. If \code{density} is \code{TRUE}, the plot will 
#' consist of multiple parts and is harder to modify.
#'
#' @examples
#' data(example_set)
#' plot_pca(example_set, color = "Injection_order", shape = "Group")
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca <- function(object, pcs = c(1, 2), all_features = FALSE, 
                     center = TRUE, scale = "uv", color = NULL,
                     shape = color, label = NULL, density = FALSE, 
                     title = "PCA", subtitle = NULL, color_scale = NA,
                     shape_scale = getOption("notame.shape_scale"), 
                     fill_scale = getOption("notame.fill_scale_dis"),
                     text_base_size = 14, point_size = 2,
                     assay.type = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(color, shape), 
                          assay.type = from)
  
  pca_results <- .pca_helper(object, pcs, center, scale, assay.type = from, ...)
  pca_scores <- pca_results$pca_scores
  pc_names <- colnames(pca_scores)
  pca_scores[color] <- colData(object)[, color]
  pca_scores[shape] <- colData(object)[, shape]
  pca_scores[label] <- colData(object)[, label]

  .scatter_plot(pca_scores, x = pc_names[1], y = pc_names[2], 
                xlab = pca_results$labels[1], ylab = pca_results$labels[2],
                color = color, shape = shape, label = label, density = density,
                title = title, subtitle = subtitle, color_scale = color_scale,
                shape_scale = shape_scale, fill_scale = fill_scale,
                text_base_size = text_base_size, point_size = point_size)
}

#' t-SNE scatter plot
#'
#' Computes t-SNE into two dimensions and plots the map points.
#' In case there are missing values, PCA is performed using the nipals method 
#' of \code{\link[pcaMethods]{pca}}, the method can be changed to "ppca" if 
#' nipals fails.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in \code{\link[pcaMethods]{prep}}. Default is 
#' '"uv" for unit variance
#' @param perplexity the perplexity used in t-SNE
#' @param pca_method the method used in PCA if there are missing values
#' @param color character, name of the column used for coloring the points. Set 
#' to NULL for black color.
#' @param shape character, name of the column used for shape. Set to NULL for 
#' uniform round shapes.
#' @param point_size numeric, size of the points.
#' @param label character, name of the column used for point labels
#' @param density logical, whether to include density plots to both axes.
#' The density curves will be split and colored by the 'color' variable.
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function.
#' Set to NA to choose the appropriate scale based on the class of the coloring 
#' variable.
#' @param shape_scale the shape scale as returned by a ggplot function
#' @param fill_scale the fill scale used for density curves.
#' If a continuous variable is used as color, density curve will be colorless.
#' @param text_base_size numeric, base size for text
#' @param point_size numeric, size of the points
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional arguments passed to \code{\link[Rtsne]{Rtsne}}
#'
#' @return A ggplot object. If \code{density} is \code{TRUE}, the plot will 
#' consist of multiple parts and is harder to modify.
#'
#' @examples
#' data(example_set)
#' plot_tsne(example_set, color = "Time", shape = "Group", perplexity = 10)
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne <- function(object, all_features = FALSE, center = TRUE, 
                      scale = "uv", perplexity = 30, pca_method = "nipals",
                      color = NULL, shape = color, label = NULL,
                      density = FALSE, title = "t-SNE",
                      subtitle = paste("Perplexity:", perplexity), 
                      color_scale = NA,
                      shape_scale = getOption("notame.shape_scale"), 
                      fill_scale = getOption("notame.fill_scale_dis"),
                      text_base_size = 14, point_size = 2, 
                      assay.type = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(color, shape),
                          assay.type = from)

  # t-SNE
  tsne_scores <- .t_sne_helper(object, center, scale, perplexity, 
                               pca_method, from, ...)
  # Add columns for plotting
  tsne_scores[color] <- colData(object)[, color]
  tsne_scores[shape] <- colData(object)[, shape]
  tsne_scores[label] <- colData(object)[, label]

  .scatter_plot(tsne_scores, x = "X1", y = "X2", color = color, shape = shape,
                label = label, density = density, title = title, 
                subtitle = subtitle, color_scale = color_scale, 
                shape_scale = shape_scale, fill_scale = fill_scale,
                text_base_size = text_base_size, point_size = point_size)
}

.mod_scatter <- function(p, data, x, y, color, label, density,
                        fill_scale, label_text_size) {
  # Add point labels
  if (!is.null(label)) {
    p <- p + ggrepel::geom_text_repel(mapping = aes(label = .data[[label]]),
                                      size = label_text_size) +
      # Remove "a" from the legend (ggrepel adds it by default)
      guides(color = guide_legend(override.aes = aes(label = ""))) 
  }

  # Add density plots to top and right
  if (density) {
    xdens <- cowplot::axis_canvas(p, axis = "x") +
      geom_density(data = data, aes(x = .data[[x]], fill = .data[[color]]),
                   alpha = 0.7, size = 0.2) +
      fill_scale

    ydens <- cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE) +
      geom_density(data = data, aes(x = .data[[y]], fill = .data[[color]]),
                   alpha = 0.7, size = 0.2) +
      coord_flip() +
      fill_scale

    p <- cowplot::insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), 
                                    position = "top")
    p <- cowplot::insert_yaxis_grob(p, ydens, grid::unit(.2, "null"), 
                                    position = "right")
    p <- cowplot::ggdraw(p)
  }
  p
}

.scatter_plot <- function(data, x, y, color, shape, label = NULL, 
                         density = FALSE, fixed = TRUE, color_scale = NA,
                         shape_scale = NULL, fill_scale = NA, title = NULL,
                         subtitle = NULL, xlab = x, ylab = y, 
                         color_lab = color, shape_lab = shape, 
                         apply_theme_bw = TRUE, 
                         text_base_size = text_base_size, 
                         point_size = point_size, 
                         label_text_size = label_text_size) {
  # Set right color scale
  if (!is.null(color_scale)) {
    if (!is.ggproto(color_scale)) {
      color_scale <- if (is.na(color_scale)) {
        if (class(data[, color]) %in% c("numeric", "integer")) {
          getOption("notame.color_scale_con")
        } else {
          getOption("notame.color_scale_dis")
        }
      } else {
        color_scale()
      }
    }
  }

  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]],
                        color = if (!is.null(color)) .data[[color]])) +
    color_scale +
    labs(title = title, subtitle = subtitle, x = xlab, 
         y = ylab, color = color_lab)
  if (apply_theme_bw) {
    p <- p + theme_bw(base_size = text_base_size)
  }
  if (fixed) {
    p <- p + coord_fixed() + theme(aspect.ratio = 1)
  }
  if (is(data[, shape], "character")) {
    data[shape] <- as.factor(data[, shape])
    warning(paste("Shape variable not given as a factor,", 
                  "converted to factor with levels",
                  paste(levels(data[, shape]), collapse = ", ")), 
            call. = FALSE)
  }
  if (is(data[, shape], "factor")) {
    if (length(levels(data[, shape])) <= 8) {
      p <- p +
        geom_point(aes(shape = .data[[shape]]), size = point_size) +
        shape_scale +
        labs(shape = shape_lab)
    } else if (is.null(shape_scale)) {
      message("Only 8 distinct shapes currently available!")
      p <- p + geom_point(size = point_size)
    }
  } else {
    p <- p + geom_point(size = point_size)
  }

  p <- .mod_scatter(p, data, x, y, color, label, density,
                   fill_scale, label_text_size)
  p
}

#' PCA loadings plot
#'
#' Computes PCA using one of the methods provided in the Bioconductor package
#' pcaMethods and plots the loadings of first principal components.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in \code{\link[pcaMethods]{prep}}. Default is 
#' "uv" for unit variance
#' @param n_features numeric vector of length two, number of top feature to plot
#' for each principal component
#' @param title,subtitle the titles of the plot
#' @param text_base_size numeric, base size for text
#' @param point_size numeric, size of the points
#' @param label_text_size numeric, size of the labels
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional arguments passed to \code{\link[pcaMethods]{prep}}
#'
#' @return A ggplot object.
#'
#' @examples
#' data(example_set)
#' plot_pca_loadings(example_set, n_features = c(2, 4))
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca_loadings <- function(object, pcs = c(1, 2), all_features = FALSE,
                              center = TRUE, scale = "uv", 
                              n_features = c(10, 10),
                              title = "PCA loadings", subtitle = NULL,
                              text_base_size = 14, point_size = 2,
                              label_text_size = 4, assay.type = NULL, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \"pcaMethods\" needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  
  .add_citation("PCA was performed using pcaMethods package:",
                citation("pcaMethods"))
  
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from, feature_ID = TRUE)

  pca_res <- pcaMethods::pca(t(assay(object, from)), nPcs = max(pcs), 
                             center = center, scale = scale, ...)

  loads <- as.data.frame(pcaMethods::loadings(pca_res))[, pcs]
  pc_names <- colnames(loads)
  loads$Feature_ID <- rownames(loads)

  features_pc1 <- loads$Feature_ID[order(abs(loads[, pc_names[1]]), 
    decreasing = TRUE)][seq_len(n_features[1])]
    
  features_pc2 <- loads$Feature_ID[order(abs(loads[, pc_names[2]]), 
    decreasing = TRUE)][seq_len(n_features[2])]

  loads <- loads[union(features_pc1, features_pc2), ]

  ggplot(loads, aes(x = .data[[pc_names[1]]], y = .data[[pc_names[2]]], 
                    label = .data[["Feature_ID"]])) +
    geom_point(size = point_size) +
    ggrepel::geom_text_repel(size = label_text_size) +
    theme_bw(base_size = text_base_size) +
    labs(title = title, subtitle = subtitle)
}


# --------------- HEXBIN PLOTS --------------------


#' PCA hexbin plot
#'
#' Computes PCA using one of the methods provided in the Bioconductor package
#' pcaMethods and plots the two first principal components as hexagonal bins, 
#' where the value of the coloring variable is summarised for each bin, by 
#' default as the mean of the values inside the bin.
#'
#' @param object a SummarizedExperiment object
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in \code{\link[pcaMethods]{prep}}. Default is 
#' "uv" for unit variance
#' @param fill character, name of the column used for coloring the hexagons
#' @param summary_fun the function used to compute the value for each hexagon
#' @param bins the number of bins in x and y axes
#' @param title,subtitle the titles of the plot
#' @param fill_scale the fill scale as returned by a ggplot function
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional arguments passed to \code{\link[pcaMethods]{pca}}
#'
#' @return A ggplot object.
#'
#' @examples
#' data(example_set)
#' plot_pca_hexbin(example_set)
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca_hexbin <- function(object, pcs = c(1, 2), all_features = FALSE, 
                            center = TRUE, scale = "uv",
                            fill = "Injection_order", summary_fun = "mean",
                            bins = 10, title = "PCA", subtitle = NULL,
                            fill_scale = getOption("notame.fill_scale_con"),
                            assay.type = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_nums = fill, assay.type = from)

  pca_results <- .pca_helper(object, pcs, center, scale, assay.type = from, ...)
  pca_scores <- pca_results$pca_scores
  pc_names <- colnames(pca_scores)
  pca_scores[fill] <- colData(object)[, fill]

  .hexbin_plot(data = pca_scores, x = pc_names[1], y = pc_names[2], 
               xlab = pca_results$labels[1], ylab = pca_results$labels[2],
               fill = fill, summary_fun = summary_fun, bins = bins, 
               fill_scale = fill_scale, title = title, subtitle = subtitle)
}

#' t-SNE hexbin plot
#'
#' Computes t-SNE into two dimensions and plots the map as hexagonal bins, 
#' where the value of the coloring variable is summarised for each bin, by 
#' default as the mean of the values inside the bin.
#' In case there are missing values, PCA is performed using the nipals method 
#' of \code{\link[pcaMethods]{pca}}, the method can be changed to "ppca" if 
#' nipals fails.
#'
#' @param object a SummarizedExperiment object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in  \code{\link[pcaMethods]{prep}}. Default is 
#' "uv" for unit variance
#' @param pca_method the method used in PCA if there are missing values
#' @param perplexity the perplexity used in t-SNE
#' @param fill character, name of the column used for coloring the hexagons
#' @param summary_fun the function used to compute the value for each hexagon
#' @param bins the number of bins in x and y axes
#' @param title,subtitle the titles of the plot
#' @param fill_scale the fill scale as returned by a ggplot function
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional arguments passed to \code{\link[Rtsne]{Rtsne}}
#'
#' @return
#' A ggplot object.
#'
#' @examples
#' data(example_set)
#' plot_tsne_hexbin(example_set, perplexity = 10)
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne_hexbin <- function(object, all_features = FALSE, center = TRUE, 
                             scale = "uv", pca_method = "nipals", 
                             perplexity = 30, fill = "Injection_order",
                             summary_fun = "mean", bins = 10, title = "t-SNE",
                             subtitle = paste("Perplexity:", perplexity),
                             fill_scale = getOption("notame.fill_scale_con"),
                             assay.type = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_nums = fill, assay.type = from)
  # t-SNE
  tsne_scores <- .t_sne_helper(object, center, scale,
                               perplexity, pca_method, from, ...)
  # Add columns for plotting
  tsne_scores[fill] <- colData(object)[, fill]

  .hexbin_plot(tsne_scores, x = "X1", y = "X2", fill = fill, 
               summary_fun = summary_fun, bins = bins,
               fill_scale = fill_scale,
               title = title, subtitle = subtitle)
}


.hexbin_plot <- function(data, x, y, fill, summary_fun = "mean", bins = 10,
                        fill_scale = NULL, title = NULL, subtitle = NULL, 
                        xlab = x, ylab = y, fill_lab = fill) {

  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], z = .data[[fill]])) +
    stat_summary_hex(bins = bins, fun = summary_fun) +
    theme_bw() +
    fill_scale +
    labs(title = title, subtitle = subtitle, x = xlab, 
         y = ylab, fill = fill_lab) +
    coord_fixed() +
    theme(aspect.ratio = 1)
  p
}


# --------- ARROW PLOTS -----------

.arrow_plot <- function(data, x, y, color, time, subject, alpha, arrow_style,
                       color_scale, title, subtitle, xlab, ylab,
                       text_base_size, line_width) {
  data <- data[order(data[, subject], data[, time]), ]

  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]], 
                        group = .data[[subject]])) +
    geom_path(arrow = arrow_style, alpha = alpha, linewidth = line_width) +
    theme_bw(base_size = text_base_size) +
    color_scale +
    coord_fixed() +
    theme(aspect.ratio = 1) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle)
  p
}

#' PCA plot with arrows
#'
#' Plots changes in PCA space according to time. All the observations of a 
#' single subject are connected by an arrow ending at the last observation.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in \code{\link[pcaMethods]{prep}}. Default is 
#' "uv" for unit variance
#' @param color character, name of the column used for coloring the arrows
#' @param time character, name of the column containing timepoints
#' @param subject character, name of the column containing subject identifiers
#' @param alpha numeric, value for the alpha parameter of the arrows 
#' (transparency)
#' @param arrow_style a description of arrow heads, the size and angle can be 
#' modified, see \code{?arrow}
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function
#' @param text_base_size the base size of the text
#' @param line_width the width of the arrows
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional arguments passed to \code{\link[pcaMethods]{pca}}
#'
#' @return A ggplot object.
#'
#' @examples
#' data(example_set)
#' plot_pca_arrows(drop_qcs(example_set), color = "Group", time = "Time",
#'   subject = "Subject_ID")
#' # If the sample size is large, plot groups separately
#' plot_pca_arrows(drop_qcs(example_set), color = "Group", 
#'                 time = "Time", subject = "Subject_ID") +
#'   facet_wrap(~Group)
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca_arrows <- function(object, pcs = c(1, 2), all_features = FALSE, 
                            center = TRUE, scale = "uv",
                            color, time, subject, alpha = 0.6,
                            arrow_style = arrow(), title = "PCA changes",
                            subtitle = NULL, 
                            color_scale = getOption("notame.color_scale_dis"),
                            text_base_size = 14, line_width = 0.5, 
                            assay.type = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(color, time, subject),
                         assay.type = from)

  pca_results <- .pca_helper(object, pcs, center, scale, assay.type = from, ...)
  pca_scores <- pca_results$pca_scores
  pc_names <- colnames(pca_scores)
  pca_scores[color] <- colData(object)[, color]
  pca_scores[time] <- colData(object)[, time]
  pca_scores[subject] <- colData(object)[, subject]

  .arrow_plot(data = pca_scores, x = pc_names[1], y = pc_names[2], 
              color = color, time = time, subject = subject,
              alpha = alpha, arrow_style = arrow_style, 
              color_scale = color_scale, title = title, subtitle = subtitle,
              xlab = pca_results$labels[1], ylab = pca_results$labels[2],
              text_base_size = text_base_size, line_width = line_width)
}

#' t-SNE plot with arrows
#'
#' Computes t-SNE into two dimensions and plots changes according to time.
#' All the observations of a single subject are connected by an arrow ending at 
#' the last observation. In case there are missing values, PCA is performed 
#' using the nipals method of \code{\link[pcaMethods]{pca}}, the method can be 
#' changed to "ppca" if nipals fails.
#'
#' @param object a SummarizedExperiment object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in  \code{\link[pcaMethods]{prep}}. Default is 
#' "uv" for unit variance
#' @param perplexity the perplexity used in t-SNE
#' @param pca_method the method used in PCA if there are missing values
#' @param color character, name of the column used for coloring the points
#' @param time character, name of the column containing timepoints
#' @param subject character, name of the column containing subject identifiers
#' @param alpha numeric, value for the alpha parameter of the arrows 
#' (transparency)
#' @param arrow_style a description of arrow heads, the size and angle can be 
#' modified, see \code{?arrow}
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function
#' @param text_base_size the base size of the text
#' @param line_width the width of the arrows
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional arguments passed to \code{\link[Rtsne]{Rtsne}}
#'
#' @return A ggplot object. If \code{density} is \code{TRUE}, the plot will 
#' consist of multiple parts and is harder to modify.
#'
#' @examples
#' data(example_set)
#' plot_tsne_arrows(drop_qcs(example_set), perplexity = 10, color = "Group", 
#'   time = "Time", subject = "Subject_ID")
#' # If the sample size is large, plot groups separately
#' plot_tsne_arrows(drop_qcs(example_set), perplexity = 10, color = "Group", 
#'   time = "Time", subject = "Subject_ID") +
#'     facet_wrap(~Group)
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne_arrows <- function(object, all_features = FALSE, center = TRUE, 
                             scale = "uv", perplexity = 30, 
                             pca_method = "nipals", color, time, subject,
                             alpha = 0.6, arrow_style = arrow(), 
                             title = "t-SNE changes",
                             subtitle = paste("Perplexity:", perplexity),
                             color_scale = getOption("notame.color_scale_dis"),
                             text_base_size = 14, line_width = 0.5,
                             assay.type = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(color, time, subject),
                         assay.type = from)

  tsne_scores <- .t_sne_helper(object, center, scale, 
                               perplexity, pca_method, from, ...)
  tsne_scores[color] <- colData(object)[, color]
  tsne_scores[time] <- colData(object)[, time]
  tsne_scores[subject] <- colData(object)[, subject]

  .arrow_plot(data = tsne_scores, x = "X1", y = "X2", color = color,
              time = time, subject = subject, alpha = alpha, 
              arrow_style = arrow_style, color_scale = color_scale,
              title = title, subtitle = subtitle, xlab = "X1", ylab = "X2",
              text_base_size = text_base_size, line_width = line_width)
}