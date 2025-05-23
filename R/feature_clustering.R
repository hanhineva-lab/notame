#' Cluster correlated features originating from the same metabolite
#'
#' Clusters features potentially originating from the same compound. Features 
#' with high Pearson correlation coefficient and small retention time 
#' difference are linked together. Then clusters are formed by setting a 
#' threshold for the relative degree that each node in a cluster needs
#' to fulfil. Each cluster is named after the feature with the highest median 
#' peak area (median abundance). This is a wrapper around numerous functions 
#' that are based on the MATLAB code by David Broadhurst.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param mz_col the column name in feature data that holds mass-to-charge 
#' ratios
#' @param rt_col the column name in feature data that holds retention times
#' @param all_features logical, should all features be included in the 
#' clustering? If FALSE, as the default, flagged features are not included in 
#' clustering
#' @param rt_window the retention time window for potential links
#' NOTE: use the same unit as the retention time
#' @param corr_thresh the correlation threshold required for potential links 
#' between features
#' @param d_thresh the threshold for the relative degree required by each node
#' @param plotting should plots be drawn for each cluster?
#' @param min_size_plotting the minimum number of features a cluster needs to 
#' have to be plotted
#' @param prefix the prefix to the files to be plotted
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return a SummarizedExperiment or MetaboSet object, with median peak area 
#' (MPA), the cluster ID, the features in the cluster, and cluster size added 
#' to feature data.
#'
#' @examples
#' data(example_set)
#' # The parameters are really weird because example data is imaginary
#' clustered <- cluster_features(example_set, rt_window = 1, corr_thresh = 0.5, 
#'   d_thresh = 0.6)
#'
#' @export
cluster_features <- function(object, mz_col = NULL, rt_col = NULL,
                             all_features = FALSE, rt_window = 1 / 60,
                             corr_thresh = 0.9, d_thresh = 0.8, 
                             plotting = TRUE, min_size_plotting = 3, 
                             prefix = NULL, assay.type = NULL) {
  # Drop flagged compounds before clustering
  from <- .get_from_name(object, assay.type)
  orig <- .check_object(object)
  object <- drop_flagged(object, all_features)
  object <- .check_object(object, check_limits = TRUE, assay.type = from,
                         feature_cols = c(mz_col, rt_col))

  if (is.null(mz_col) || is.null(rt_col)) {
    cols <- .find_mz_rt_cols(rowData(object))
  }
  mz_col <- mz_col %||% cols$mz_col
  rt_col <- rt_col %||% cols$rt_col

  data <- as.data.frame(t(assay(object, from)))
  features <- as.data.frame(rowData(object))
  # Start log
  log_text(paste("\nStarting feature clustering at", Sys.time()))

  # Find connections between features in each Split
  conn <- data.frame()
  for (s in unique(features$Split)) {
    log_text(paste("Finding connections between features in", s))
    features_tmp <- features[features$Split == s, ]
    data_tmp <- data[, features_tmp$Feature_ID]

    conn_tmp <- find_connections(data = data_tmp, features = features_tmp,
                                 corr_thresh = corr_thresh,
                                 rt_window = rt_window, name_col = "Feature_ID",
                                 mz_col = mz_col, rt_col = rt_col)
    
    conn <- rbind(conn, conn_tmp)
    log_text(paste("Found", nrow(conn_tmp), "connections in", s))
  }
  log_text(paste("Found", nrow(conn), "connections"))


  # Form clusters
  clusters <- find_clusters(conn, d_thresh)
  lens <- vapply(clusters, function(x) length(x$features), integer(1))
  log_text(paste("Found", sum(lens > 1),
                 "clusters of 2 or more features, clustering finished at",
                 Sys.time()))
  # Compute median peak area and assing cluster ID
  features$MPA <- apply(assay(object, from), 1, finite_median)
  features <- assign_cluster_id(data, clusters, features, "Feature_ID")

  if (plotting) {
    visualize_clusters(data = data, features = features, clusters = clusters,
                       min_size = min_size_plotting, rt_window = rt_window,
                       name_col = "Feature_ID", mz_col = mz_col, 
                       rt_col = rt_col, file_path = prefix)
    log_text(paste("Saved cluster plots to:", prefix))
  }
  # Add cluster IDs to the ORIGINAL object (flagged features still there)
  clustered <- join_rowData(orig, features[c("Feature_ID", "MPA",
                                            "Cluster_ID", "Cluster_size",
                                            "Cluster_features")])
  if (!is.null(attr(clustered, "original_class"))) {
    clustered <- as(clustered, "MetaboSet")
    attr(object, "original_class") <- NULL
  }                                          
  clustered
}

#' Assign Cluster ID to features
#'
#' Assigns a cluster ID to all features that are part of a cluster with 2 or 
#' more features.
#'
#' @param data data frame of the original LC-MS data
#' @param clusters a list of clusters as returned by find_clusters
#' @param features data frame with feature information, feature data
#' @param name_col character, name of the column in features that contains 
#' feature names
#'
#' @inherit find_connections return examples
#'
#' @return A data frame similar to features, with cluster ID added.
#'
#' @noRd
assign_cluster_id <- function(data, clusters, features, name_col) {
  if (!"MPA" %in% colnames(features)) {
    features$MPA <- vapply(data[, features[, name_col]], 
                           finite_median, numeric(1))
  }

  features$Cluster_ID <- features[, name_col]
  features$Cluster_features <- features[, name_col]
  features$Cluster_size <- 1
  n_clust <- length(clusters)
  for (cluster in clusters) {
    if (length(cluster$features) > 1) {
      # Which features are in the cluster
      idx <- features[, name_col] %in% cluster$features
      # The cluster is named for the feature with the largest median peak area
      features_tmp <- features[idx, ]
      max_mpa_idx <- which(features_tmp$MPA == max(features_tmp$MPA, 
                                                   na.rm = TRUE))[1]
      max_mpa_feature <- features_tmp[max_mpa_idx, name_col]
      # Saving some information about the clusters
      features$Cluster_ID[idx] <- paste0("Cluster_", max_mpa_feature)
      features$Cluster_size[idx] <- length(cluster$features)
      features$Cluster_features[idx] <- paste(sort(cluster$features), 
                                              collapse = ";")
    }
  }
  features
}

#' Compress clusters of features to a single feature
#'
#' This function compresses clusters found by cluster_features, keeping only 
#' the feature with the highest median peak area. The features that were 
#' discarded are recorded in feature data, under Cluster_features.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#'
#' @return A SummarizedExperiment or MetaboSet object with only one feature per 
#' cluster.
#'
#' @examples
#' data(example_set)
#' clustered <- cluster_features(example_set, 
#'   rt_window = 1, corr_thresh = 0.5, d_thresh = 0.6)
#' compressed <- compress_clusters(clustered)
#'
#' @seealso \code{\link{cluster_features}}
#'
#' @export
compress_clusters <- function(object) {
  object <- .check_object(object)
  cluster_names <- rowData(object)$Cluster_ID
  if (is.null(cluster_names)) {
    stop("No 'Cluster_ID' found in rowData(object), ",
         "please run cluster_features first!")
  }
  # Get only "real" clusters
  clusters <- cluster_names[grepl("^Cluster_", cluster_names)] %>%
    gsub("^Cluster_", "", .) %>%
    unique()
  alone_features <- cluster_names[!grepl("^Cluster_", cluster_names)]
  # This ensures the order of the features stays the same
  idx <- rowData(object)$Feature_ID %in% c(clusters, alone_features)

  object <- object[idx, ]
  log_text(paste("Clusters compressed, left with", nrow(object), "features"))
  if (!is.null(attr(object, "original_class"))) {
    object <- as(object, "MetaboSet")
    attr(object, "original_class") <- NULL
  }
  object
}

#' Extract information of the features in clusters
#'
#' For each cluster, the LC-MS data of the feature with largest median peak 
#' area is retained, all the features inside every cluster are recorded.
#'
#' @param data data frame of the original LC-MS data
#' @param features data frame holding the feature information
#' @param name_col name_col character, name of the column in features that
#' contains feature names
#'
#' @inherit find_connections return examples
#'
#' @return A list of two items:
#' \itemize{
#' \item cdata: a new data frame with the combined LC-MS data
#' \item cfeatures: data frame, feature information per cluster
#' }
#'
#' @noRd
pull_clusters <- function(data, features, name_col) {
  cluster_names <- features$Cluster_ID
  if (is.null(cluster_names)) {
    stop("No 'Cluster_ID' found in features, ",
         "please run assign_cluster_id first!")
  }
  # Get only "real" clusters
  clusters <- cluster_names[grepl("^Cluster_", cluster_names)] %>%
    gsub("^Cluster_", "", .) %>%
    unique()
  alone_features <- cluster_names[!grepl("^Cluster_", cluster_names)]
  # This ensures the order of the features stays the same
  idx <- features[, name_col] %in% c(clusters, alone_features)
  # Get sample features
  sample_cols <- setdiff(colnames(data), features[, name_col])

  features <- features[idx, ]
  data <- data[, c(sample_cols, features[, name_col])]

  return(list(cdata = data, cfeatures = features))
}


#' Find out which features are correlated within a specified retention time 
#' window
#'
#' A part of the peak clustering algorithm. Iterates over all possible pairs of 
#' features and records a connection between them if a) they have a Pearson 
#' correlation coefficient higher than \code{corr_thresh} and b) their 
#' retention time difference is less than \code{rt_window}. 
#'
#' @param data data frame with the abundances of features, with features as 
#' columns
#' @param features data frame with feature information, feature data
#' @param corr_thresh numeric, the threshold of correlation to use in linking 
#' features
#' @param rt_window numeric, the retention time window to use in linking 
#' features. NOTE: you need to use the same unit as in the retention time column
#' @param name_col character, name of the column in features that contains 
#' feature names
#' @param mz_col character, name of the column in features that contains mass-
#' to-charge ratios
#' @param rt_col character, name of the column in features that contains 
#' retention times
#'
#' @examples 
#' \dontshow{.old_wd <- setwd(tempdir())}
#' data(example_set)
#' data <- combined_data(example_set)
#' features <- rowData(example_set)
#' features$MPA <- sapply(data[, features[, "Feature_ID"]], finite_median)
#' conn <- find_connections(data = data, features = features,
#'   corr_thresh = 0.4, rt_window = 2,
#'   name_col = "Feature_ID", mz_col = "Mass", rt_col = "RetentionTime")
#' clusters <- find_clusters(connections = conn, d_thresh = 0.6)
#' features_clustered <- assign_cluster_id(data, clusters, features,
#'   name_col = "Feature_ID")
#' visualize_clusters(data, features, clusters, min_size = 3,
#'   rt_window = 2, name_col = "Feature_ID", mz_col = "Mass",
#'   rt_col = "RetentionTime", file_path = "./clusters")
#' pulled <- pull_clusters(data, features_clustered, name_col = "Feature_ID")
#' \dontshow{setwd(.old_wd)}
#'
#' @return A data frame of pairs of signals that are linked together:
#' \itemize{
#' \item x & y: indexes and names of the signals
#' \item cor: correlation coefficient
#' \item mz_diff & rt_diff: mass and retention time difference
#' }
#'
#' @noRd
find_connections <- function(data, features, corr_thresh = 0.9,
                             rt_window = 1 / 60, name_col, mz_col, rt_col) {
  d <- data[features[, name_col]]
  if (ncol(d) < 2) {
    stop("Need at least 2 features to do any clustering!")
  }
  n <- nrow(features)
  
  connections <- BiocParallel::bplapply(seq_len(n -1), function(i) {
    if (i %% 100 == 0) {
      message(i)
    }
    connections_tmp <- data.frame()
    for (j in (i + 1):n) {
      rt_diff <- features[j, rt_col] - features[i, rt_col]
      cor_coef <- stats::cor(d[, i], d[, j], use = "na.or.complete")
      if (!is.na(cor_coef)) {
        if (abs(rt_diff) < rt_window && cor_coef > corr_thresh) {
          mz_diff <- features[j, mz_col] - features[i, mz_col]
          connections_tmp <- rbind(connections_tmp, 
                                   data.frame(x = features[i, name_col], 
                                              y = features[j, name_col],
                                              cor = cor_coef, 
                                              rt_diff = rt_diff, 
                                              mz_diff = mz_diff))
        }
      }
    }
    connections_tmp
  })
  connections <- do.call(rbind, connections)
}

#' Extract the densely connected clusters
#'
#' First forms clusters of compounds that are linked together. Then the 
#' clusters are pruned so that in the final clusters, each feature is linked to 
#' at least a set percentage of the other features in the cluster.
#'
#' @param connections data frame of pairs of signals that are linked together,
#' output of find_connections
#' @param d_thresh numeric, the minimum degree required for each signal in a 
#' cluster expressed as a percentage of the maximum degree in the cluster
#'
#' @inherit find_connections return examples
#'
#' @return  A list of clusters, each a list of:
#' \itemize{
#' \item features: character vector of the names of the features included in 
#' the cluster
#' \item graph: an igraph object of the cluster
#' }
#'
#' @noRd
find_clusters <- function(connections, d_thresh = 0.8) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package \"igraph\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation(paste0("igraph package was used to construct networks of", 
                       " features for feature clustering:"), 
                citation("igraph"))

  # Construct graph from the given edges
  g <- igraph::graph_from_edgelist(as.matrix(connections[seq_len(2)]), 
                                   directed = FALSE)
  g <- igraph::set.edge.attribute(graph = g, name = "weight", 
                                  value = connections$cor)

  # Initialize list of clusters
  clusters <- list()
  k <- 1

  # Repeatedly extract densely connected clusters from the graph
  while (length(igraph::V(g))) {
    # Connected components of the remaining graph
    comp <- igraph::decompose(g)
    n_comp <- length(comp)
    message(n_comp, " components found")

    # Only keep the densely connected part of each component (subgraph)
    clusters_tmp <- BiocParallel::bplapply(comp, function(subg) {
      n_nodes <- length(igraph::V(subg))
      d <- igraph::degree(subg)
      # The limit of the degree a node needs to be kept
      d_lim <- round(d_thresh * (n_nodes - 1))
      
      if (n_nodes >= 3) {
        # Remove the node with the smallest degree until all nodes in the
        # cluster have a degree above the limit
        while (any(d < d_lim)) {
          idx <- which(d == min(d))
          if (length(idx) > 1) {
            edgesums <- igraph::strength(subg, vids = igraph::V(subg)$name[idx])
            idx <- idx[which(edgesums == min(edgesums))[1]]
          }
          subg <- igraph::delete.vertices(subg, v = igraph::V(subg)[idx])
          d <- igraph::degree(subg)
          n_nodes <- n_nodes - 1
          d_lim <- round(d_thresh * (n_nodes - 1))
        }
      }

      # Record the final cluster and remove the nodes from the main graph
      list(list(features = names(igraph::V(subg)),
                graph = subg))
    })
    clusters_tmp <- do.call(c, clusters_tmp)

    for (j in seq_along(clusters_tmp)) {
      clusters[[k]] <- clusters_tmp[[j]]
      k <- k + 1
      subg <- clusters_tmp[[j]]$graph
      g <- igraph::delete.vertices(g, v = names(igraph::V(subg)))
    }
  }
  clusters
}


# A helper function for plotting, scales the values in X
# between new min and max
.rescale <- function(x, new_min, new_max) {
  y <- (new_max - new_min) * (x - min(x)) / (max(x) - min(x)) + new_min
  # If all MPAs are equal, the sizes are NaN (possibly other reasons)
  if (sum(is.na(y))) {
    y <- rep(mean(c(new_min, new_max)), length(x))
  }
  y
}

# UNFINISHED!!
.plot_graph <- function(features, cluster, name_col, mz_col, rt_col) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package \'igraph\' needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  
  .add_citation(paste0("igraph package was used to construct networks of ",
                       "features for feature clustering:"), 
                citation("igraph"))

  # Ensure a correct order of the rows
  g <- cluster$graph
  vertices <- data.frame(Name = igraph::V(g)$name, stringsAsFactors = FALSE)
  colnames(vertices) <- name_col
  features_tmp <- dplyr::left_join(vertices, features, by = name_col)

  # Scaling of MPA to correct size
  # Square root to scale area, not radius
  size <- sqrt(.rescale(features_tmp$MPA, new_min = 15^2, new_max = 40^2))

  if (length(igraph::E(g)) <= 20) {
    edge_labels <- as.character(round(igraph::E(g)$weight, digits = 2))
  } else {
    edge_labels <- NULL
  }

  g$palette <- RColorBrewer::brewer.pal(
    n = max(3, length(unique(igraph::degree(g)))), 
    name = "Blues")
  plot(g, vertex.label = as.character(features_tmp[, mz_col]),
       vertex.size = size, vertex.label.dist = 0.1 * size,
       vertex.label.degree = -pi / 2, vertex.color = igraph::degree(g),
       edge.label = edge_labels)
}

.plot_features <- function(features, cluster, name_col, 
                          mz_col, rt_col, rt_window) {

  features_tmp <- features[features[, name_col] %in% cluster$features, ]
  p1 <- ggplot(features_tmp, aes(.data[[mz_col]], .data[["MPA"]])) +
    geom_point(size = 3, color = "steelblue4") +
    geom_segment(aes(x = .data[[mz_col]], yend = .data[["MPA"]], 
                     xend = .data[[mz_col]]),
                 y = 0, color = "steelblue4") +
    ggrepel::geom_label_repel(aes(label = .data[[mz_col]]), 
                              color = "steelblue4") +
    theme_minimal() +
    xlim(0.9 * min(features_tmp[, mz_col], na.rm = TRUE),
         1.15 * max(features_tmp[, mz_col], na.rm = FALSE)) +
    expand_limits(y = 0) +
    labs(x = "Mass-to-charge ratio", y = "Median Peak Area")

  features_tmp$rtmin <- features_tmp[, rt_col] - rt_window
  features_tmp$rtmax <- features_tmp[, rt_col] + rt_window

  p2 <- ggplot(features_tmp, aes(.data[[rt_col]], .data[[mz_col]])) +
    geom_point(size = 3, color = "steelblue4") +
    geom_errorbarh(aes(xmin = .data$rtmin, xmax = .data$rtmax), 
                   color = "steelblue4") +
    theme_minimal() +
    labs(x = "Retention time", y = "Mass-to-charge ratio", 
         title = "Retention time & tolerance")

  plot(p1)
  plot(p2)
}

.plot_heatmaps <- function(data, features, cluster, name_col, mz_col, rt_col) {
  features_tmp <- features[features[, name_col] %in% cluster$features, ]

  n <- length(cluster$features)
  mz_rt <- data.frame()

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mz_rt <- rbind(mz_rt, data.frame(
        x = features_tmp[i, name_col],
        y = features_tmp[j, name_col],
        mz_diff = features_tmp[i, mz_col] - features_tmp[j, mz_col],
        rt_diff = features_tmp[i, rt_col] - features_tmp[j, rt_col],
        stringsAsFactors = FALSE))
    }
  }

  mz_ord <- features_tmp[, name_col][order(features_tmp[, mz_col])]
  mz_rt$x <- factor(mz_rt$x, levels = mz_ord)
  mz_rt$y <- factor(mz_rt$y, levels = rev(mz_ord))

  p1 <- ggplot(mz_rt, aes(x = .data$x, y = .data$y, fill = .data$mz_diff)) +
    geom_tile(color = "grey80") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
    scale_fill_gradient2()

  if (nrow(mz_rt) <= 10) {
    p1 <- p1 + geom_text(aes(label = round(.data$mz_diff, digits = 2)))
  }

  plot(p1)
}


#' Visualize clusters of features
#'
#' Draws multiple visualizations of each cluster, creating a separate file for
#' each cluster.
#'
#' @param data data frame with the abundances of features, with features as 
#' columns
#' @param features data frame with feature information, feature data
#' @param clusters a list of clusters as returned by find_clusters
#' @param min_size the minimum number of features a cluster needs to have to be 
#' plotted
#' @param rt_window numeric, the retention time window to use in linking 
#' features. NOTE you need to use the same unit as in the retention time column
#' @param name_col character, name of the column in features that contains 
#' feature names
#' @param mz_col character, name of the column in features that contains
#' mass-to-charge ratios
#' @param rt_col character, name of the column in features that contains 
#' retention times
#' @param file_path the prefix to the files to be plotted
#'
#' @inherit find_connections return examples
#'
#' @noRd
visualize_clusters <- function(data, features, clusters, min_size, rt_window,
                               name_col, mz_col, rt_col, file_path) {
  for (i in seq_along(clusters)) {
    if (i %% 100 == 0) {
      message(i, " / ", length(clusters))
    }
    cluster <- clusters[[i]]
    if (length(cluster$features) >= min_size) {
      features_tmp <- features[features[, name_col] %in% cluster$features, ]
      cluster_id <- features_tmp$Cluster_ID[1]
      grDevices::pdf(paste0(file_path, cluster_id, ".pdf"), 
                     width = 10, height = 10)
      .plot_heatmaps(data, features, cluster, name_col, mz_col, rt_col)
      .plot_features(features, cluster, name_col, mz_col, rt_col, rt_window)
      .plot_graph(features, cluster, name_col, mz_col, rt_col)
      grDevices::dev.off()
    }
  }
}
