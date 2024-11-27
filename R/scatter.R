# ------- HELPER FUNCTIONS ----------------

# Helper function for computing PCA
.pca_helper <- function(object, pcs, center, scale, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \'pcaMethods\' needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("PCA was performed using pcaMethods package:",
                citation("pcaMethods"))
  res_pca <- pcaMethods::pca(t(assay(object)), nPcs = max(pcs), scale = scale, 
                             center = center, ...)
  pca_scores <- as.data.frame(pcaMethods::scores(res_pca))[, pcs]
  r2 <- summary(res_pca)["R2", pcs]
  labels <- paste0(paste0("PC", pcs), " (", scales::percent(r2), ")")

  return(list(pca_scores = pca_scores, labels = labels))
}

# Helper function for computing t-SNE
.t_sne_helper <- function(object, center, scale, perplexity, pca_method, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \'pcaMethods\' needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Package \'Rtsne\' needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  .add_citation("Rtsne package was used for t-SNE figures:", citation("Rtsne"))
  prepd <- pcaMethods::prep(t(assay(object)), center = center, scale = scale)

  if (sum(is.na(prepd)) > 0) {
    res_pca <- pcaMethods::pca(t(assay(object)), method = pca_method, 
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
#' @param object a SummarizedExperiment or MetaboSet object
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit 
#' variance
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
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return A ggplot object. If \code{density} is \code{TRUE}, the plot will 
#' consist of multiple parts and is harder to modify.
#'
#' @examples
#' plot_pca(example_set, color = "Injection_order", shape = "Group")
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca <- function(object, pcs = c(1, 2), all_features = FALSE, 
                     center = TRUE, scale = "uv", color = group_col(object),
                     shape = color, label = NULL, density = FALSE, 
                     title = "PCA", subtitle = NULL, color_scale = NA,
                     shape_scale = getOption("notame.shape_scale"), 
                     fill_scale = getOption("notame.fill_scale_dis"),
                     text_base_size = 14, point_size = 2, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  object <- check_object(object)

  pca_results <- .pca_helper(object, pcs, center, scale, ...)
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
#' of \code{pcaMethods::pca}, the  method can be changed to "ppca" if nipals 
#' fails.
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit 
#' variance
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
#' @param ... additional arguments passed to \code{Rtsne::Rtsne}
#'
#' @return A ggplot object. If \code{density} is \code{TRUE}, the plot will 
#' consist of multiple parts and is harder to modify.
#'
#' @examples
#' plot_tsne(example_set, color = "Time", shape = "Group", perplexity = 10)
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne <- function(object, all_features = FALSE, center = TRUE, 
                      scale = "uv", perplexity = 30, pca_method = "nipals",
                      color = group_col(object), shape = color, label = NULL,
                      density = FALSE, title = "t-SNE",
                      subtitle = paste("Perplexity:", perplexity), 
                      color_scale = NA,
                      shape_scale = getOption("notame.shape_scale"), 
                      fill_scale = getOption("notame.fill_scale_dis"),
                      text_base_size = 14, point_size = 2, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  object <- check_object(object)

  # t-SNE
  tsne_scores <- .t_sne_helper(object, center, scale, 
                               perplexity, pca_method, ...)
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
#' @param object a SummarizedExperiment or MetaboSet object
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit 
#' variance
#' @param n_features numeric vector of length two, number of top feature to plot
#' for each principal component
#' @param title,subtitle the titles of the plot
#' @param text_base_size numeric, base size for text
#' @param point_size numeric, size of the points
#' @param label_text_size numeric, size of the labels
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return A ggplot object.
#'
#' @examples
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
                              label_text_size = 4, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \"pcaMethods\" needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  
  .add_citation("PCA was performed using pcaMethods package:",
                citation("pcaMethods"))
  
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  object <- check_object(object)

  pca_res <- pcaMethods::pca(t(assay(object)), nPcs = max(pcs), 
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
#' @param object a SummarizedExperiment or MetaboSet object
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit 
#' variance
#' @param fill character, name of the column used for coloring the hexagons
#' @param summary_fun the function used to compute the value for each hexagon
#' @param bins the number of bins in x and y axes
#' @param title,subtitle the titles of the plot
#' @param fill_scale the fill scale as returned by a ggplot function
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return A ggplot object.
#'
#' @examples
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
                            ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  object <- check_object(object)

  pca_results <- .pca_helper(object, pcs, center, scale, ...)
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
#' of \code{pcaMethods::pca}, the  method can be changed to "ppca" if nipals 
#' fails.
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit 
#' variance
#' @param pca_method the method used in PCA if there are missing values
#' @param perplexity the perplexity used in t-SNE
#' @param fill character, name of the column used for coloring the hexagons
#' @param summary_fun the function used to compute the value for each hexagon
#' @param bins the number of bins in x and y axes
#' @param title,subtitle the titles of the plot
#' @param fill_scale the fill scale as returned by a ggplot function
#' @param ... additional arguments passed to Rtsne::Rtsne
#'
#' @return
#' A ggplot object.
#'
#' @examples
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
                             ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  object <- check_object(object)

  # t-SNE
  tsne_scores <- .t_sne_helper(object, center, scale,
                               perplexity, pca_method, ...)
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
#' @param object a SummarizedExperiment or MetaboSet object
#' @param pcs numeric vector of length 2, the principal components to plot
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit 
#' variance
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
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return A ggplot object.
#'
#' @examples
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
                            color = group_col(object), time = time_col(object),
                            subject = subject_col(object), alpha = 0.6,
                            arrow_style = arrow(), title = "PCA changes",
                            subtitle = NULL, 
                            color_scale = getOption("notame.color_scale_dis"),
                            text_base_size = 14, line_width = 0.5, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  object <- check_object(object)

  pca_results <- .pca_helper(object, pcs, center, scale, ...)
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
#' using the nipals method of \code{pcaMethods::pca}, the method can be changed 
#' to "ppca" if nipals fails.
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually 
#' yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit 
#' variance
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
#' @param ... additional arguments passed to \code{Rtsne::Rtsne}
#'
#' @return A ggplot object. If \code{density} is \code{TRUE}, the plot will 
#' consist of multiple parts and is harder to modify.
#'
#' @examples
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
                             pca_method = "nipals",
                             color = NULL,
                             time = NULL, 
                             subject = NULL,
                             alpha = 0.6, arrow_style = arrow(), 
                             title = "t-SNE changes",
                             subtitle = paste("Perplexity:", perplexity),
                             color_scale = getOption("notame.color_scale_dis"),
                             text_base_size = 14, line_width = 0.5, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  object <- check_object(object)

  tsne_scores <- .t_sne_helper(object, center, scale, 
                               perplexity, pca_method, ...)
  tsne_scores[color] <- colData(object)[, color]
  tsne_scores[time] <- colData(object)[, time]
  tsne_scores[subject] <- colData(object)[, subject]

  .arrow_plot(data = tsne_scores, x = "X1", y = "X2", color = color,
              time = time, subject = subject, alpha = alpha, 
              arrow_style = arrow_style, color_scale = color_scale,
              title = title, subtitle = subtitle, xlab = "X1", ylab = "X2",
              text_base_size = text_base_size, line_width = line_width)
}


# -------- VOLCANO PLOT -----------

minus_log10 <- scales::trans_new("minus_log10", 
                                 transform = function(x) -log10(x),
                                 inverse = function(x) 10 ^ -x)

#' Volcano plot
#'
#' Draws a volcano plot of effect size and p-values.
#'
#' @param object a SummarizedExperiment object or a data frame. If x is a 
#' SummarizedExperiment object, feature data of x is used. If x is a data 
#' frame, it is used as is.
#' @param x,p the column names of effect size (x-axis) and p-values
#' @param p_fdr column name of FDR corrected p-values, used to draw a line 
#' showing the fdr-corrected significance level
#' @param color column name used to color the plots
#' @param p_breaks a numerical vector of the p_values to show on the y-axis
#' @param fdr_limit the significance level used in the experiment
#' @param log2_x logical, whether effect size should be plotted on a log2 axis.
#' @param center_x_axis logical, whether x-axis should be centered. If 
#' \code{TRUE}, the "zero-effect" will be on the middle of the plot. The 
#' "zero effect" is 0 if \code{log2_x = FALSE} and 1 if \code{log2_x = TRUE}
#' @param x_lim numerical vector of length 2 for manually setting the x-axis 
#' limits
#' @param label column name used to label the plots
#' @param label_limit numeric, p-value which is used to limit label plotting. 
#' Defaults to 0.05.
#' @param color_scale the color scale as returned by a ggplot function
#' @param title,subtitle the title and subtitle of the plot
#' @param text_base_size integer, base size for text in figures
#' @param label_text_size numeric, size of the labels
#' @param ...  parameters passed to \code{\link[ggplot2]{geom_point}},
#' such as shape and alpha values. New aesthetics can
#' also be passed using \code{mapping = aes(...)}.
#'
#' @return A ggplot object.
#'
#' @examples
#' # naturally, this looks messy as there are not enough p-values
#' lm_results <- perform_lm(drop_qcs(example_set), 
#'   formula_char = "Feature ~ Group")
#' volcano_plot(lm_results,
#'   x = "GroupB_Estimate",
#'   p = "GroupB_P", p_fdr = "GroupB_P_FDR",
#'   label = "Feature_ID",
#'   fdr_limit = 0.1
#' )
#'
#' @export
setGeneric("volcano_plot", signature = "object",
  function(object, x, p, p_fdr = NULL, color = NULL,
          p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
          log2_x = FALSE, center_x_axis = TRUE, x_lim = NULL, 
          label = NULL, label_limit = 0.05, 
          color_scale = getOption("notame.color_scale_con"),
          title = "Volcano plot", subtitle = NULL,
          text_base_size = 14, label_text_size = 4, ...) {
    standardGeneric("volcano_plot")
  }
)

#' @rdname volcano_plot
#' @export
setMethod("volcano_plot", c(object = "SummarizedExperiment"),
  function(object, x, p, p_fdr = NULL, color = NULL,
           p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
           log2_x = FALSE, center_x_axis = TRUE, x_lim = NULL, 
           label = NULL, label_limit = 0.05, 
           color_scale = getOption("notame.color_scale_con"),
           title = "Volcano plot", subtitle = NULL,
           text_base_size = 14, label_text_size = 4, ...) {
    .volcano_plotter(rowData(object), x, p, p_fdr, color, p_breaks, fdr_limit,
                     log2_x, center_x_axis, x_lim, label, label_limit,
                     color_scale, title, subtitle,
                     text_base_size, label_text_size, ...)
  }
)

#' @rdname volcano_plot
#' @export
setMethod("volcano_plot", c(object = "data.frame"),
  function(object, x, p, p_fdr = NULL, color = NULL,
           p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
           log2_x = FALSE, center_x_axis = TRUE, x_lim = NULL, label = NULL,
           label_limit = 0.05,
           color_scale = getOption("notame.color_scale_con"),
           title = "Volcano plot", subtitle = NULL,
           text_base_size = 14, label_text_size = 4, ...) {
    .volcano_plotter(object, x, p, p_fdr, color, p_breaks, fdr_limit,
                     log2_x, center_x_axis, x_lim, label, label_limit,
                     color_scale, title, subtitle, text_base_size,
                     label_text_size, ...)
  }
)


.volcano_plotter <- function(data, x, p, p_fdr, color, p_breaks, fdr_limit, 
                            log2_x, center_x_axis, x_lim, label, label_limit,
                            color_scale, title, subtitle, text_base_size,
                            label_text_size, ...) {
                              
                              
  if (center_x_axis && !is.null(x_lim)) {
    warning("Manually setting x-axis limits overrides x-axis centering.")
    center_x_axis <- FALSE
  }
  if (min(data[, p]) > max(p_breaks)) {
    warning("All the p-values are larger than the p-value breaks supplied.", 
            " Consider using larger p_breaks for plotting.")
  }

  pl <- ggplot(data, aes(x = .data[[x]], y = .data[[p]],
                         color = if (!is.null(color)) .data[[color]])) +
    geom_point(...) +
    color_scale +
    theme_bw(base_size = text_base_size) +
    labs(title = title, subtitle = subtitle, color = color) +
    # Only show the specified p-values, which might be unevenly spaced
    theme(panel.grid.minor.y = element_blank(), 
          axis.ticks.y = element_blank())

  if (!is.null(p_fdr)) {
    if (any(data[, p_fdr] < fdr_limit)) {
      # Add horizontal line with the FDR < 0.05 limit
      q_limit <- max(data[data[, p_fdr] < fdr_limit, p], na.rm = TRUE)
      # sec_axis writes e.g. "q < 0.05" on the right side of the plot
      pl <- pl +
        geom_hline(yintercept = q_limit, linetype = "dashed") +
        scale_y_continuous(trans = minus_log10, breaks = p_breaks,
                           labels = as.character(p_breaks),
                           sec.axis = sec_axis(~., breaks = q_limit, 
                                               labels = paste("q =",fdr_limit)))
    } else {
      warning("None of the FDR-adjusted p-values are below the significance", 
              " level, not plotting the horizontal line.", call. = FALSE)
      pl <- pl +
        scale_y_continuous(trans = minus_log10, breaks = p_breaks,
                           labels = as.character(p_breaks))
    }
  } else {
    pl <- pl +
      scale_y_continuous(trans = minus_log10, breaks = p_breaks,
                         labels = as.character(p_breaks))
  }

  if (log2_x) {
    if (center_x_axis) {
      x_lim <- max(abs(log2(data[, x])))
      x_lim <- c(2^(-x_lim), 2^x_lim)
    }
    pl <- pl +
      scale_x_continuous(trans = "log2", limits = x_lim)
  } else {
    if (center_x_axis) {
      x_lim <- max(abs(data[, x]))
      x_lim <- c(-x_lim, x_lim)
    }
    pl <- pl +
      scale_x_continuous(limits = x_lim)
  }

  if (!is.null(label)) {
    if (label %in% colnames(data)) {
      label_data <- data[data[, p] < label_limit, ]
      pl <- pl +
        ggrepel::geom_label_repel(data = label_data, 
                                  mapping = aes(label = .data[[label]]),
                                  seed = 313, alpha = 0.5, 
                                  size = label_text_size,
                                  force = 10, max.overlaps = 50)
    } else {
      warning("Label column not found, not plotting them")
    }
  }

  pl
}


# ---------- MANHATTAN PLOT -------


#' Manhattan plot
#'
#' Draws a (directed) Manhattan plot of p-values and versus e.g. retention time 
#' or mass-to-charge ratio. If effect size and direction is supplied, the -
#' log10(p-value) on the y-axis will be multiplied by the direction (sign) of 
#' the effect, so part of the points will "drop" from the p = 1 (-log10(p) = 0) 
#' line. This results in a so-called directed Manhattan plot.
#'
#' @param object a SummarizedExperiment object or a data frame. If x is a 
#' SummarizedExperiment object, feature data of x is used. If x is a data 
#' frame, it is used as is.
#' @param x,p the column names of x-axis and p-values
#' @param effect column name of effect size (should have negative and positive 
#' values).
#' @param p_fdr column name of FDR corrected p-values, used to draw a line 
#' showing the fdr-corrected significance level
#' @param color column name used to color the plots
#' @param p_breaks a numerical vector of the p_values to show on the y-axis
#' @param fdr_limit the significance level used in the experiment
#' @param x_lim,y_lim numerical vectors of length 2 for manually setting the 
#' axis limits
#' @param color_scale the color scale as returned by a ggplot function
#' @param title,subtitle the title and subtitle of the plot
#' @param ...  parameters passed to \code{\link[ggplot2]{geom_point}},
#' such as shape and alpha values. New aesthetics can
#' also be passed using \code{mapping = aes(...)}.
#'
#' @return A ggplot object.
#'
#' @examples
#' # naturally, this looks messy as there are not enough p-values
#' lm_results <- perform_lm(drop_qcs(example_set), 
#'   formula_char = "Feature ~ Group")
#' lm_data <- dplyr::left_join(as.data.frame(rowData(example_set)), lm_results)
#' # Traditional Manhattan plot from data frame
#' manhattan_plot(lm_data,
#'   x = "Average_Mz",
#'   p = "GroupB_P", p_fdr = "GroupB_P_FDR",
#'   fdr_limit = 0.1
#' )
#' # Directed Manhattan plot from SummarizedExperiment
#' with_results <- join_rowData(example_set, lm_results)
#' manhattan_plot(with_results,
#'   x = "Average_Mz", effect = "GroupB_Estimate",
#'   p = "GroupB_P", p_fdr = "GroupB_P_FDR",
#'   fdr_limit = 0.1
#' )
#'
#' @export
setGeneric("manhattan_plot", signature = "object",
  function(object, x, p, effect = NULL, p_fdr = NULL, color = NULL,
           p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
           x_lim = NULL, y_lim = NULL,
           color_scale = getOption("notame.color_scale_con"),
           title = "Manhattan plot", subtitle = NULL, ...) {
    standardGeneric("manhattan_plot")
  }
)

#' @rdname manhattan_plot
#' @export
setMethod("manhattan_plot", c(object = "SummarizedExperiment"),
  function(object, x, p, effect = NULL, p_fdr = NULL, color = NULL,
           p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
           x_lim = NULL, y_lim = NULL,
           color_scale = getOption("notame.color_scale_con"),
           title = "Manhattan plot", subtitle = NULL, ...) {
    .manhattan_plotter(as.data.frame(rowData(object)), x, p, effect, p_fdr, 
                       color, p_breaks, fdr_limit, x_lim, y_lim, color_scale, 
                       title, subtitle, ...)
  }
)

#' @rdname manhattan_plot
#' @export
setMethod("manhattan_plot", c(object = "data.frame"),
  function(object, x, p, effect = NULL, p_fdr = NULL, color = NULL,
           p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
           x_lim = NULL, y_lim = NULL,
           color_scale = getOption("notame.color_scale_con"),
           title = "Manhattan plot", subtitle = NULL, ...) {
    .manhattan_plotter(object, x, p, effect, p_fdr, color, p_breaks, fdr_limit,
                       x_lim, y_lim, color_scale, title, subtitle, ...)
  }
)


.manhattan_plotter <- function(data, x, p, effect, p_fdr, color, p_breaks,
                              fdr_limit, x_lim, y_lim, color_scale,
                              title, subtitle, ...) {
  if (min(data[, p]) > max(p_breaks)) {
    warning("All the p-values are larger than the p-value breaks supplied.", 
            " Consider using larger p_breaks for plotting.")
  }

  if (!is.null(effect)) {
    data$y <- -log10(data[, p]) * sign(data[, effect])
    p_labels <- outer(c(-1, 1), p_breaks) %>%
      as.vector() %>%
      as.character()
    p_breaks <- outer(c(-1, 1), -log10(p_breaks)) %>% as.vector()
    p_labels <- p_labels[order(p_breaks)]
    p_breaks <- sort(p_breaks)
  } else {
    data$y <- -log10(data[, p])
    p_labels <- as.character(p_breaks)
    p_breaks <- -log10(p_breaks)
    p_labels <- p_labels[order(p_breaks)]
    p_breaks <- sort(p_breaks)
  }

  pl <- ggplot(data, aes(x = .data[[x]], y = .data[["y"]],
                         color = if (!is.null(color)) .data[[color]])) +
    geom_point(...) +
    color_scale +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          axis.ticks.y = element_blank()) +
    geom_hline(yintercept = 0, color = "grey") +
    labs(title = title, subtitle = subtitle, y = "p-value", color = color)


  if (!is.null(p_fdr)) {
    if (any(data[, p_fdr] < fdr_limit)) {
      # Add horizontal line with the FDR < 0.05 limit
      q_limit <- max(data[data[, p_fdr] < fdr_limit, p], na.rm = TRUE)
      # sec_axis writes e.g. "q < 0.05" on the right sifde of the plot
      if (!is.null(effect)) {
        pl <- pl +
          geom_hline(yintercept = log10(q_limit), linetype = "dashed") +
          geom_hline(yintercept = -log10(q_limit), linetype = "dashed") +
          scale_y_continuous(
            breaks = p_breaks, labels = p_labels, limits = y_lim,
            sec.axis = sec_axis(~., breaks = c(log10(q_limit), -log10(q_limit)),
                                labels = rep(paste("q =", fdr_limit), 2)))
      } else {
        pl <- pl +
          geom_hline(yintercept = -log10(q_limit), linetype = "dashed") +
          scale_y_continuous(
            breaks = p_breaks, labels = p_labels, limits = y_lim,
            sec.axis = sec_axis(~., breaks = -log10(q_limit), 
                                labels = paste("q =", fdr_limit)))
      }
    } else {
      warning("None of the FDR-adjusted p-values are below the significance", 
              " level, not plotting the horizontal line.", call. = FALSE)
      pl <- pl +
        scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim)
    }
  } else {
    pl <- pl +
      scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim)
  }

  pl
}


# -----M/Z vs RT PLOT -------

#' Plot m/z vs retention time plot (cloud plot)
#'
#' Plots a scatter plot of results of statistical tests, where each point 
#' represents a feature. The plot has retention time on x-axis, m/z on y-axis 
#' and the size of the points is scaled based on p-value
#'
#' @param object a SummarizedExperiment object or a data frame. If x is a 
#' SummarizedExperiment object, rowData(x) is used.
#' If x is a data frame, it is used as is.
#' @param p_col the column name containing p-values. This is used to scale the 
#' size of the points.
#' @param p_limit numeric, limits plotted features by p-values. If NULL, plots 
#' all features.
#' @param mz_col,rt_col the column names for m/z and retention time. If NULL, 
#' automatic detection is attempted.
#' @param color the column name used to color the points
#' @param title The plot title
#' @param subtitle The plot subtitle
#' @param color_scale color scale as returned by a ggplot function. Defaults to 
#' current continuous color scale.
#' @param all_features logical, should all features be retained? Should be used 
#' only if x is a SummarizedExperiment object.
#' @param ...  parameters passed to \code{\link[ggplot2]{geom_point}},
#' such as shape and alpha values. New aesthetics can
#' also be passed using \code{mapping = aes(...)}.
#'
#' @return A ggplot object.
#'
#' @examples
#' # Compute results from a linear model
#' lm_results <- perform_lm(example_set, formula_char = "Feature ~ Group")
#' with_results <- join_rowData(example_set, lm_results)
#'
#' # Plot from the SummarizedExperiment object
#' # automatically facet by analytical mode in variable Split
#' mz_rt_plot(with_results, p_col = "GroupB_P", color = "GroupB_Estimate")
#'
#' # Plot the results from the results dataframe
#' mz_rt_plot(with_results, p_col = "GroupB_P", color = "GroupB_Estimate")
#'
#' @export
setGeneric("mz_rt_plot", signature = "object",
  function(object, p_col = NULL, p_limit = NULL, mz_col = NULL, rt_col = NULL,
           color = NULL, title = "m/z retention time", subtitle = NULL,
           color_scale = getOption("notame.color_scale_con"), ...) {
    standardGeneric("mz_rt_plot")
  }
)

#' @rdname mz_rt_plot
#' @export
setMethod("mz_rt_plot", c(object = "SummarizedExperiment"),
  function(object, p_col = NULL, p_limit = NULL, mz_col = NULL, rt_col = NULL,
           color = NULL, title = "m/z vs retention time", subtitle = NULL,
           color_scale = getOption("notame.color_scale_con"), 
           all_features = FALSE) {
    .mz_rt_plotter(as.data.frame(rowData(drop_flagged(object, all_features))),
                   p_col, p_limit, mz_col, rt_col, color, title, subtitle, 
                   color_scale, all_features)
  }
)

#' @rdname mz_rt_plot
#' @export
setMethod("mz_rt_plot", c(object = "data.frame"),
  function(object, p_col = NULL, p_limit = NULL, mz_col = NULL, rt_col = NULL,
           color = NULL, title = "m/z vs retention time", subtitle = NULL,
           color_scale = getOption("notame.color_scale_con")) {
    .mz_rt_plotter(object, p_col, p_limit, mz_col, rt_col, color, 
                   title, subtitle, color_scale)
  }
)


.mz_rt_plotter <- function(x, p_col, p_limit, mz_col, rt_col, color, 
                          title, subtitle, color_scale, all_features) {
  if (!is.null(p_limit) && !is.null(p_col)) {
    x <- x[which(x[, p_col] < p_limit), ]
    if (nrow(x) == 0) {
      stop("No features with p-values smaller than ", p_limit, " found.")
    }
    message("All features with p-values larger than ", p_limit, " dropped.\n")
  }

  if (is.null(mz_col) || is.null(rt_col)) {
    mz_rt_cols <- .find_mz_rt_cols(x)
    mz_col <- mz_col %||% mz_rt_cols$mz_col
    rt_col <- rt_col %||% mz_rt_cols$rt_col
  }

  p <- ggplot(x, aes(x = .data[[rt_col]], y = .data[[mz_col]], 
                     size = .data[[p_col]], 
                     color = if (!is.null(color)) .data[[color]])) +
    geom_point(alpha = 0.6) +
    scale_size_continuous(trans = minus_log10, range = c(0.5, 3),
                          breaks = c(0.05, 0.01, 0.001, 1e-4), 
                          labels = as.character(c(0.05, 0.01, 0.001, 1e-4))) +
    theme_bw() +
    color_scale +
    labs(title = title, subtitle = subtitle, color = color, 
         x = "Retention time", y = "Mass-to-charge ratio", size = "p-value") +
    # Scales for m/z and rt
    scale_x_continuous(breaks = seq(0, ceiling(max(x[, rt_col])))) +
    scale_y_continuous(breaks = seq(0, ceiling(max(x[, mz_col])), 250))

  # If multiple splits are given, plot them separately
  if (length(unique(x$Split)) > 1) {
    message("Multiple splits detected, plotting them to separate panes.\n")
    p <- p +
      facet_wrap(~Split, dir = "v")
  }

  p
}
