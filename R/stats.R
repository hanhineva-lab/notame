## Helper function for calculating summary statistics
.calc_summary <- function(f_levels, grouping_cols, groups, group_names, 
                          funs, result_row_template) {
  result_row <- result_row_template
  for (fname in names(funs)) {
    tmp <- tapply(f_levels, groups, funs[[fname]])
    if (is.null(grouping_cols[1])) {
      result_row[fname] <- tmp[1]
    } else {
      result_row[paste(group_names, fname, sep = "_")] <- tmp
    }
  }
  return(result_row)
}

#' Summary statistics
#'
#' Computes summary statistics for each feature, possibly grouped by a factor.
#' The statistics include mean, standard deviation (sd), median,
#' median absolute deviation (mad), minimum (min), maximum (max)
#' as well as 25% and 75% quantiles (Q25 & Q75).
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param grouping_cols character vector, the columns by which grouping should 
#' be done. Use \code{NA} to compute statistics without grouping.
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @examples
#' data(example_set)
#' # Group by "Group"
#' sum_stats <- summary_statistics(example_set, grouping_cols = "Group")
#' # Group by Group and Time
#' sum_stats <- summary_statistics(example_set, 
#'   grouping_cols = c("Group", "Time"))
#' # No Grouping
#' sum_stats <- summary_statistics(example_set)
#'
#' @return A data frame with the summary statistics.
#'
#' @export
summary_statistics <- function(object, grouping_cols = NULL,
                               assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(grouping_cols), 
                         assay.type = from)
  data <- combined_data(object, from)
  features <- rownames(object)
  # Get sample grouping and group names for saving results accordingly
  if (is.null(grouping_cols)[1]) {
    groups <- rep(1, nrow(data))
    group_names <- ""
  } else {
    # Single grouping column
    if (length(grouping_cols) == 1) {
      groups <- data[, grouping_cols]
      if (is(groups, "factor")) {
        group_names <- levels(groups)
      } else {
        group_names <- unique(groups)
      }
    } else {
      groups <- rep("", nrow(data))
      for (grouping_col in grouping_cols) {
        tmp_group <- paste(grouping_col, data[, grouping_col], sep = "_")
        groups <- paste(groups, tmp_group, sep = "_")
      }
      groups <- as.factor(gsub("^_", "", groups))
      group_names <- levels(groups)
    }
  }
  # Define functions to use
  funs <- list(mean = finite_mean, sd = finite_sd, median = finite_median,
               mad = finite_mad, min = finite_min,
               Q25 = function(x) finite_quantile(x, probs = 0.25),
               Q75 = function(x) finite_quantile(x, probs = 0.75),
               max = finite_max)
  # Initialize named vector for the results
  result_row_template <- rep(0, times = length(group_names) * length(funs))
  if (!is.null(grouping_cols[1])) {
    var_names <- expand.grid(names(funs), group_names)
    names(result_row_template) <- paste(var_names$Var2, var_names$Var1, 
                                        sep = "_")
  } else {
    names(result_row_template) <- names(funs)
  }
  # Calculate statistics
  statistics <- t(data.frame(
    BiocParallel::bplapply(data[, features], .calc_summary, 
                           grouping_cols = grouping_cols, groups = groups,
                           group_names = group_names, funs = funs,
                           result_row_template = result_row_template), 
    stringsAsFactors = FALSE))
  # Convert to data frame
  data.frame(Feature_ID = rownames(statistics), statistics, row.names = NULL)
}

#' Statistics cleaning
#'
#' Uses regexp to remove unnecessary columns from statistics results data frame.
#' Can also rename columns effectively.
#'
#' @param df data frame, statistics results
#' @param remove list, should contain strings that are matching to unwanted 
#' columns
#' @param rename named list, names should contain matches that are replaced 
#' with values
#' @param summary logical, should summary columns be added
#' @param p_limit numeric, limit for p-values to be counted
#' @param fdr logical, should summary be done with fdr-fixed values
#'
#' @return A data frame with removed and/or renamed columns.
#'
#' @examples
#' data(example_set)
#' # Simple manipulation to linear model results
#' lm_results <- perform_lm(drop_qcs(example_set), 
#'   formula_char = "Feature ~ Group + Time")
#' lm_results <-summarize_results(lm_results,
#'   rename = c("GroupB" = "GroupB_vs_A", "Time2" = "Time2_vs_1")
#' )
#'
#' @export
summarize_results <- function(df, remove = c("Intercept", "CI95", "Std_error",
                                             "t_value", "z_value", "R2"),
                              rename = NULL, summary = TRUE, p_limit = 0.05,
                              fdr = TRUE) {
  df <- df[, !grepl(paste(remove, collapse = "|"), colnames(df))]
  if (!is.null(rename)) {
    for (name in names(rename)) {
      colnames(df) <- gsub(name, rename[name], colnames(df))
    }
  }
  if (summary) {
    ifelse(fdr, p_cols <- colnames(df)[grep("_P_FDR$", colnames(df))],
           p_cols <- colnames(df)[grep("_P$", colnames(df))])
    df$Low_p_values <- apply(df[, p_cols], 1, 
                             function(x) sum(x < p_limit, na.rm = TRUE))
    df$Lowest_p <- apply(df[, p_cols], 1, 
                         function(x) finite_min(x))
    df$Column <- apply(df[, p_cols], 1, 
      function(x) {
        m <- finite_min(x)
        p_cols[which(x == m)]
      })
  }
  df
}

.calc_cohens_d <- function(feature, group1, group2) {
  f1 <- group1[, feature]
  f2 <- group2[, feature]
  d <- data.frame(
    Feature_ID = feature,
    Cohen_d = (finite_mean(f2) - finite_mean(f1)) /
      sqrt((finite_sd(f1)^2 + finite_sd(f2)^2) / 2),
    stringsAsFactors = FALSE)
}

.help_cohens_d <- function(object, group, id, time, assay.type) {
  data <- combined_data(object, assay.type)
  features <- rownames(object)
  group_levels <- levels(data[, group])
  time_levels <- NULL

  if (is.null(time)) {
    group1 <- data[which(data[, group] == group_levels[1]), ]
    group2 <- data[which(data[, group] == group_levels[2]), ]
    log_text(paste("Starting to compute Cohen's D between groups",
                   paste(rev(group_levels), collapse = " & ")))
  } else {
    time_levels <- levels(data[, time])
    # Split to time points
    time1 <- data[which(data[, time] == time_levels[1]), ]
    time2 <- data[which(data[, time] == time_levels[2]), ]
    common_ids <- intersect(time1[, id], time2[, id])
    rownames(time1) <- time1[, id]
    rownames(time2) <- time2[, id]
    time1 <- time1[common_ids, ]
    time2 <- time2[common_ids, ]
    if (!identical(time1[, group], time2[, group])) {
      stop("Groups of subjects do not match between time points.", 
           call. = FALSE)
    }
    # Change between time points
    new_data <- time2[, features] - time1[, features]
    # Split to groups
    group1 <- new_data[which(time1[, group] == levels(time1[, group])[1]), ]
    group2 <- new_data[which(time1[, group] == levels(time1[, group])[2]), ]

    log_text(paste("Starting to compute Cohen's D between groups",
                   paste(rev(group_levels), collapse = " & "),
                   "from time change",
                   paste(rev(time_levels), collapse = " - ")))
  }
  
  ds <- BiocParallel::bplapply(features, FUN = .calc_cohens_d,
                               group1 = group1, group2 = group2)
  
  ds <-  do.call(rbind, ds)

  if (is.null(time_levels)) {
    colnames(ds)[2] <- 
      paste0(group_levels[2], "_vs_", group_levels[1], "_Cohen_d")
  } else {
    colnames(ds)[2] <- paste0(group_levels[2], "_vs_", group_levels[1],
                              "_", time_levels[2], "_minus_", 
                              time_levels[1], "_Cohen_d")
  }

  log_text("Cohen's D computed.")
  ds
}


.help_cohens_d_time <- function(res, group_combos, object, group,
                                id, time, assay.type) {
  count_obs_geq_than <- function(x, n) {
    sum(x >= n)
  }
  if (is.null(id)) {
    stop("Please specify id column.", call. = FALSE)
  }
  time_combos <- utils::combn(levels(colData(object)[, time]), 2)
  for (i in seq_len(ncol(group_combos))) {
    for (j in seq_len(ncol(time_combos))) {
      object_split <- object[, which(
        colData(object)[, group] %in% group_combos[, i] &
          colData(object)[, time] %in% time_combos[, j])]
      colData(object_split) <- droplevels(colData(object_split))
      # Check data is valid for Cohen's D
      group_table <- table(colData(object_split)[, c(id, group)])
      time_table <- table(colData(object_split)[, c(id, time)])
      column <- paste0("Cohen_d_", group_combos[1, i], "_", group_combos[2, i],
                       "_", time_combos[2, j], "_minus_", time_combos[1, j])
      if (any(apply(group_table, 2, count_obs_geq_than, 2) < 2)) {
        warning("In ", column, ": Groups don't have two observations of",
                "at least two subjects, skipping!")
        next
      }
      if (any(apply(time_table, 1, count_obs_geq_than, 2) != 0)) {
        warning("In ", column,": Same subject recorded more than once at",
                " same time, skipping!")
        next
      }
      if (any(apply(group_table, 1, count_obs_geq_than, 1) != 1)) {
        warning("In ", column, 
                ": Same subject recorded in two groups, skipping!")
        next
      }
      if (!all(apply(time_table, 1, count_obs_geq_than, 1) != 1)) {
        warning("One or more subject(s) missing time points, ", column, 
                " will be counted using common subjects in time points!")
      }
      if (is.null(res)) {
        res <- .help_cohens_d(object_split, group, id, time, assay.type)
      } else {
        res <- dplyr::full_join(res,
          .help_cohens_d(object_split, group, id, time, assay.type),
          by = "Feature_ID")
      }  
    }
  }
  res
}


#' Cohen's D
#'
#' Computes Cohen's D for each feature. If time and ID are supplied,
#' change between two time points is computed for each subject,
#' and Cohen's d is computed from the changes.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param id character, name of the subject ID column
#' @param group character, name of the group column
#' @param time character, name of the time column
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A data frame with Cohen's d for each feature.
#'
#' @examples
#' data(example_set)
#' d_results <- cohens_d(drop_qcs(example_set), group = "Group")
#' d_results_time <- cohens_d(drop_qcs(example_set),
#'   group = "Group", time = "Time", id = "Subject_ID"
#' )
#'
#' @export
cohens_d <- function(object, group, id = NULL,
                     time = NULL, assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_factors = c(group, time),
                          pheno_chars = id, pheno_cols = id, 
                          assay.type = from, feature_ID = TRUE)
  res <- NULL

  group_combos <- utils::combn(levels(colData(object)[, group]), 2)

  if (is.null(time)) {
    for (i in seq_len(ncol(group_combos))) {
      object_split <- object[, which(
        colData(object)[, group] %in% c(group_combos[1, i], 
                                        group_combos[2, i]))]
      colData(object_split) <- droplevels(colData(object_split))

      if (is.null(res)) {
        res <- .help_cohens_d(object_split, group, id, time, from)
      } else {
        res <- dplyr::full_join(res,
          .help_cohens_d(object_split, group, id, time, from),
          by = "Feature_ID")
      }
    }
  } else {
    res <- .help_cohens_d_time(res, group_combos, object, group, id, time, from)
  }
  rownames(res) <- res$Feature_ID
  res
}

.calc_fold_change <- function(feature, group, data, groups) {
  result_row <- rep(NA_real_, ncol(groups))
  # Calculate fold changes
  tryCatch({
    for (i in seq_len(ncol(groups))) {
      group1 <- data[data[, group] == groups[1, i], feature]
      group2 <- data[data[, group] == groups[2, i], feature]
      result_row[i] <- finite_mean(group2) / finite_mean(group1)
    }
  })
  result_row
}

#' Fold change
#'
#' Computes fold change between each group for each feature.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param group character, name of the group column
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A data frame with fold changes for each feature.
#'
#' @examples
#' data(example_set)
#' # Between groups
#' fc <- fold_change(example_set, group = "Group")
#' # Between time points
#' fc <- fold_change(example_set, group = "Time")
#'
#' @export
fold_change <- function(object, group, assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = group, assay.type = from)
  log_text("Starting to compute fold changes.")
  data <- combined_data(object, from)
  groups <- utils::combn(levels(data[, group]), 2)
  features <- rownames(object)
  # Calculate fold changes between groupings
  results_df <- BiocParallel::bplapply(features, FUN = .calc_fold_change,
                                       group, data, groups)
  results_df <- do.call(rbind, results_df)
  # Create comparison labels for result column names
  comp_labels <- groups %>%
    t() %>%
    as.data.frame() %>%
    tidyr::unite("Comparison", "V2", "V1", sep = "_vs_")
  comp_labels <- paste0(comp_labels[, 1], "_FC")
  results_df <- data.frame(features, results_df, stringsAsFactors = FALSE)
  colnames(results_df) <- c("Feature_ID", comp_labels)
  rownames(results_df) <- results_df$Feature_ID

  log_text("Fold changes computed.")

  # Order the columns accordingly
  results_df[c("Feature_ID", comp_labels[order(comp_labels)])]
}


.calc_correlation <- function(var_pair, data1, data2, id, ...) {
  x_tmp <- var_pair["x", ]
  y_tmp <- var_pair["y", ]
  cor_tmp <- NULL
  
  tryCatch(
    {
      if (is.null(id)) {
        cor_tmp <- stats::cor.test(data1[, x_tmp], data2[, y_tmp], ...)
      } else {
        id_tmp <- data1[, id]
        cor_tmp <- data.frame(id_var = id_tmp, x_var = data1[, x_tmp],
                              y_var = data2[, y_tmp]) %>%
        rmcorr::rmcorr(participant = .$id_var, measure1 = .$x_var,
                       measure2 = .$y_var, dataset = .)
        cor_tmp <- list(estimate = cor_tmp$r, p.value = cor_tmp$p)
      }
    },
    error = function(e) message(x_tmp, " vs", y_tmp, ": ", e$message)
  )
  if (is.null(cor_tmp)) {
    cor_tmp <- list(estimate = NA, p.value = NA)
  }
  data.frame(X = x_tmp, Y = y_tmp, Correlation_coefficient = cor_tmp$estimate,
             Correlation_P = cor_tmp$p.value, stringsAsFactors = FALSE)
}


.help_correlation_tests <- function(var_pairs, data1, data2, 
                                    id, fdr, duplicates , ...) {
  # Prepare pairs for bplapply iteration
  var_pairs <- apply(var_pairs, 1, data.frame)
  # Calculate correlations for each pair
  cor_results <- BiocParallel::bplapply(var_pairs, .calc_correlation, 
                                        data1, data2, id, ...)
  cor_results <- do.call(rbind, cor_results)
  
  if (duplicates) {
    cor_results_dup <- cor_results
    cor_results_dup$X <- cor_results$Y
    cor_results_dup$Y <- cor_results$X
    # Remove possible duplicated correlations of a variable with itself
    cor_results_dup <- dplyr::filter(cor_results_dup, .data$X != .data$Y)
    cor_results <- rbind(cor_results, cor_results_dup)
  }

  # FDR correction
  if (fdr) {
    flags <- rep(NA_character_, nrow(cor_results))
    cor_results <- .adjust_p_values(cor_results, flags)
  }

  rownames(cor_results) <- seq_len(nrow(cor_results))
  cor_results
}

#' Correlation test
#'
#' Performs a correlation test between two sets of variables. All the variables 
#' must be either feature names or column names of pheno data (sample 
#' information).
#' There are two ways to use this function:
#' either provide a set of variables as \code{x}, and all correlations between
#' those variables are computed. Or provide two distinct sets of variables 
#' \code{x, y} and correlations between each x variable
#' and each y variable are computed.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param x character vector, names of variables to be correlated
#' @param y character vector, either identical to x (the default) or a distinct 
#' set of variables to be correlated against x
#' @param id character, column name for subject IDs. If provided, the 
#' correlation will be computed using the rmcorr package
#' @param object2 optional second object. If provided, x 
#' variables will be taken from object and y variables will be taken from 
#' object2. Both objects should have the same number of samples.
#' @param fdr logical, whether p-values from the correlation test should be 
#' adjusted with FDR correction
#' @param all_pairs logical, whether all pairs between x and y should be tested.
#' If FALSE, x and y give the exact pairs of variables to test, and should have 
#' the same length.
#' @param duplicates logical, whether correlations should be duplicated. If 
#' \code{TRUE}, each correlation will be included in the results twice, where 
#' the order of the variables '(which is x and which is y) is changed. Can be 
#' useful for e.g. plotting a heatmap of the results, see examples of
#' \code{\link{plot_effect_heatmap}}
#' @param assay.type1 character, assay of object(1) to be used in case of 
#' multiple assays
#' @param assay.type2 character, assay of object2 to be used in case of 
#' multiple assays
#' @param ... other parameters passed to \code{\link{cor.test}}, such as method
#'
#' @return A data frame with the results of correlation tests: the pair of 
#' variables, correlation coefficient and p-value.
#'
#' @examples
#' data(example_set)
#' # Correlations between all features
#' correlations <- perform_correlation_tests(example_set, 
#'   x = rownames(example_set), id = "Subject_ID")
#'
#' # Spearman Correlations between features and sample information variables
#' # Drop QCs and convert time to numeric
#' no_qc <- drop_qcs(example_set)
#' no_qc$Time <- as.numeric(no_qc$Time)
#' correlations <- perform_correlation_tests(no_qc,
#'   x = rownames(example_set),
#'   y = c("Time", "Injection_order"), method = "spearman"
#' )
#'
#' # Correlations between variables from two distinct objects
#' cross_object_cor <- perform_correlation_tests(hilic_neg_sample,
#'   x = rownames(hilic_neg_sample),
#'   object2 = hilic_pos_sample,
#'   y = rownames(hilic_pos_sample),
#'   all_pairs = FALSE
#' )
#' @seealso \code{\link{cor.test}}, \code{\link[rmcorr]{rmcorr}}
#'
#' @export
perform_correlation_tests <- function(object, x, y = x, id = NULL, 
                                      object2 = NULL, fdr = TRUE, 
                                      all_pairs = TRUE, duplicates = FALSE,
                                      assay.type1 = NULL, assay.type2 = NULL,
                                      ...) {
  log_text("Starting correlation tests.")
  
  from1 <- .get_from_name(object, assay.type1)
  object <- .check_object(object, pheno_cols = id, assay.type = from1)
  data1 <- combined_data(object, from1)

  if (!is.null(object2)) {
    if (ncol(object) != ncol(object2)) {
      stop("The objects have different numbers of samples")
    }
    from2 <- .get_from_name(object2, assay.type2)
    object2 <- .check_object(object2, pheno_factors = id, assay.type = from2)
    data2 <- combined_data(object2, from2)
    log_text("Performing correlation tests for two objects")
  } else {
    if (!is.null(assay.type2)) {
      stop("When using a single object, assay.type2 is not considered.")
    } else {
      data2 <- data1
      log_text("Performing correlation tests for single object")
    }
  }

  # Checks for repeated measures correlation
  if (!is.null(id)) {
    if (!requireNamespace("rmcorr", quietly = TRUE)) {
      stop("Package \'rmcorr\' needed for this function to work.",
           " Please install it.", call. = FALSE)
    }
    if (!id %in% colnames(data1) || !id %in% colnames(data2)) {
      stop("id column not found.", call. = FALSE)
    }
    if (!identical(data1[, id], data2[, id])) {
      stop("ids do not match between the two objects:",
           " make sure the subjects are in the same order!", call. = FALSE)
    }
    .add_citation(paste0("rmcorr package was used to compute correlations with",
                         "repeated measurements:"),
                  citation("rmcorr"))
  }

  # All x and y should be columns names of combined data
  not_found <- setdiff(x, colnames(data1))
  not_found <- c(not_found, setdiff(y, colnames(data2)))
  if (length(not_found)) {
    stop(paste("Following variables do not match to know variables in",
               "the object(s):", paste(not_found, collapse = ", ")))
  }
  
  if (all_pairs) {
    # If the same variable is present in x and y, the correlation would be 
    # computed twice. This makes sure only unique combinations are treated.
    if (identical(x, y)) {
      var_pairs <- utils::combn(x, 2) %>%
        t() %>%
        data.frame(stringsAsFactors = FALSE)
      colnames(var_pairs) <- c("x", "y")
      # Add correlations of all variables with themselves (useful for plotting)
      var_pairs <- rbind(var_pairs, data.frame(x = x, y = x, 
                                               stringsAsFactors = FALSE))
    } else if (is.null(object2) && length(intersect(x, y))) {
      stop("Currently only identical x & y or completely separate x & y are", 
           " supported for one object")
    } else {
      var_pairs <- expand.grid(x, y, stringsAsFactors = FALSE)
      colnames(var_pairs) <- c("x", "y")
    }
  } else {
    if (length(x) != length(y)) {
      stop("If all_pairs = FALSE, x and y should have the same length")
    }
    var_pairs <- data.frame(x = x, y = y, stringsAsFactors = FALSE)
  }

  cor_results <- .help_correlation_tests(var_pairs, data1, data2,
                                         id, fdr, duplicates, ...)
  log_text("Correlation tests performed.")
  cor_results
}


.calc_auc <- function(feature, data, pheno_data, time, subject, group) {
  result_row <- rep(NA_real_, nrow(pheno_data))
  # Compute AUC for each subject in each group
  tryCatch({
    for (j in seq_len(nrow(pheno_data))) {
      subset_idx <- data[, subject] == pheno_data[j, subject] & 
        data[, group] == pheno_data[j, group]
      result_row[j] <- PK::auc(time = as.numeric(data[subset_idx, time]),
                               conc = data[subset_idx, feature], 
                               design = "complete")$est[1]
    }
  })
  matrix(result_row, nrow = 1, dimnames = list(feature, pheno_data$Sample_ID))
}

#' Area under curve
#'
#' Compute area under curve (AUC) for each subject and feature.
#' Creates a pseudo SummarizedExperiment or MetaboSet object, where the 
#' "samples" are subjects (or subject/group combinations in case the same 
#' subjects are submitted to different treatments) and the "abundances" are 
#' AUCs. This object can then be  used to compute results of e.g. t-tests of 
#' AUCs between groups.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param time,subject,group column names of pheno data holding time, 
#' subject and group labels
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A pseudo SummarizedExperiment or MetaboSet object with the AUCs.
#'
#' @examples
#' data(example_set)
#' # Drop QC samples before computing AUCs
#' aucs <- perform_auc(drop_qcs(example_set), time = "Time", 
#'                     subject = "Subject_ID", group = "Group")
#' # t-test with the AUCs
#' t_test_results <- perform_t_test(aucs, formula_char = "Feature ~ Group")
#'
#' @seealso \code{\link[PK]{auc}}
#'
#' @export
perform_auc <- function(object, time, subject, group, assay.type = NULL) {
  if (!requireNamespace("PK", quietly = TRUE)) {
    stop("Package \"PK\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  .add_citation("PK package was used to compute AUC:", citation("PK"))

  log_text("Starting AUC computation")
  
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_factors = c(time, group),
                          pheno_chars = subject, assay.type = from)
  data <- combined_data(object, from)

  # Create new pheno data, only one row per subject and group
  pheno_data <- data[, c(subject, group)] %>%
    dplyr::distinct() %>%
    tidyr::unite("Sample_ID", subject, group, remove = FALSE)
  # QC and Injection_order columns to pass validObject check
  pheno_data$QC <- "Sample"
  pheno_data$Injection_order <- seq_len(nrow(pheno_data))
  rownames(pheno_data) <- pheno_data$Sample_ID

  features <- rownames(object)
  aucs <- BiocParallel::bplapply(features, .calc_auc, data, pheno_data,
                                 time, subject, group)
  aucs <- do.call(rbind, aucs)

  # Construct new SummarizedExperiment object (with all modes together)
  new_object <- SummarizedExperiment(assays = aucs, 
                                     rowData = rowData(object),
                                     colData = pheno_data) %>%
    merge_objects()

  log_text("AUC computation finished")
  
  if (!is.null(attr(object, "original_class"))) {
    new_object <- as(new_object, "MetaboSet")
  }
  new_object
}

# Helper function for FDR correction
.adjust_p_values <- function(x, flags) {
  p_cols <- colnames(x)[grepl("_P$", colnames(x))]
  for (p_col in p_cols) {
    p_values <- x[, p_col, drop = TRUE]
    p_values[!is.na(flags)] <- NA
    x <- tibble::add_column(.data = x,
                            FDR = stats::p.adjust(p_values, method = "BH"),
                            .after = p_col)
    p_idx <- which(colnames(x) == p_col)
    colnames(x)[p_idx + 1] <- paste0(p_col, "_FDR")
  }
  x
}

# Helper function for filling missing rows in results files with NAs
# Some statistical tests may fail for some features, due to e.g. missing values.
.fill_results <- function(results_df, features) {
  # Add NA rows for features where the test failed
  results_df <- results_df %>% dplyr::select("Feature_ID", dplyr::everything())
  missing_features <- setdiff(features, results_df$Feature_ID)
  fill_nas <- matrix(NA, nrow = length(missing_features), 
                     ncol = ncol(results_df) - 1) %>%
    as.data.frame()
  results_fill <- data.frame(Feature_ID = missing_features, fill_nas)
  rownames(results_fill) <- missing_features
  colnames(results_fill) <- colnames(results_df)
  results_df <- rbind(results_df, results_fill) %>% as.data.frame()
  rownames(results_df) <- results_df$Feature_ID
  # Set Feature ID to the original order
  results_df <- results_df[features, ]
  results_df
}

.help_perform_test <- function(feature, data, formula_char, result_fun, ...) {
  # Replace "Feature" with the current feature name
  tmp_formula <- gsub("Feature", feature, formula_char)
  # Run test
  result_row <- result_fun(feature = feature, 
                           formula = stats::as.formula(tmp_formula), 
                           data = data, ...)
  # In case Feature is used as predictor, make the column names match
  if (!is.null(result_row)) {
    colnames(result_row) <- gsub(feature, "Feature", colnames(result_row))
  }
  result_row
}


# Helper function for running a variety of simple statistical tests
.perform_test <- function(object, formula_char, result_fun, all_features, 
                          fdr = TRUE, packages = NULL, assay.type, ...) {
  data <- combined_data(object, assay.type)
  features <- rownames(object)

  results_df <- BiocParallel::bplapply(features, .help_perform_test, data,
                                       formula_char, result_fun, ...)
                                       
  results_df <- dplyr::bind_rows(results_df)

  if (nrow(results_df) == 0) {
    stop("All the tests failed.",
         "To see the problems, run the tests without parallelization.", 
         call. = FALSE)
  }
  # Rows full of NA for features where the test failed
  results_df <- .fill_results(results_df, features)

  # FDR correction
  if (fdr) {
    if (all_features) {
      flags <- rep(NA_character_, nrow(results_df))
    } else {
      flags <- flag(object)
    }
    results_df <- .adjust_p_values(results_df, flags)
  }
  results_df
}

#' Linear models
#'
#' Fits a linear model separately for each feature. Returns all relevant
#' statistics.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param formula_char character, the formula to be used in the linear model 
#' (see Details)
#' @param all_features should all features be included in FDR correction?
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional parameters passed to \code{\link{lm}}
#'
#' @return A data frame with one row per feature, with all the
#' relevant statistics of the linear model as columns.
#'
#' @details The linear model is fit on combined_data(object). Thus, column names
#' in pheno data can be specified. To make the formulas flexible, the word 
#' "Feature" must be used to signal the role of the features in the formula. 
#' "Feature" will be replaced by the actual Feature IDs during model fitting, 
#' see the example.
#'
#' @examples
#' data(example_set)
#' # A simple example without QC samples
#' # Features predicted by Group and Time
#' lm_results <- perform_lm(drop_qcs(example_set), 
#'   formula_char = "Feature ~ Group + Time")
#'
#' @seealso \code{\link[stats]{lm}}
#'
#' @export
perform_lm <- function(object, formula_char, all_features = FALSE, 
                       assay.type = NULL, ...) {
  log_text("Starting linear regression.")
  
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from)

  lm_fun <- function(feature, formula, data) {
    # Try to fit the linear model
    fit <- NULL
    tryCatch(
      {
        fit <- stats::lm(formula, data = data, ...)
      },
      error = function(e) message(feature, ": ", e$message))
    if (is.null(fit) || sum(!is.na(data[, feature])) < 2) {
      result_row <- NULL
    } else {
      # Gather coefficients and CIs to one data frame row
      coefs <- summary(fit)$coefficients
      confints <- stats::confint(fit, level = 0.95)
      coefs <- data.frame(Variable = rownames(coefs), coefs, 
                          stringsAsFactors = FALSE)
      confints <- data.frame(Variable = rownames(confints), confints,
                             stringsAsFactors = FALSE)

      result_row <- dplyr::left_join(coefs, confints, by = "Variable") %>%
        dplyr::rename("Std_Error" = "Std..Error", "t_value" = "t.value",
                      "P" = "Pr...t..", "LCI95" = "X2.5..", 
                      "UCI95" = "X97.5..") %>%
        tidyr::gather("Metric", "Value", -"Variable") %>%
        tidyr::unite("Column", "Variable", "Metric", sep = "_") %>%
        tidyr::spread("Column", "Value")
      # Add R2 statistics and feature ID
      result_row$R2 <- summary(fit)$r.squared
      result_row$Adj_R2 <- summary(fit)$adj.r.squared
      result_row$Feature_ID <- feature
    }
    result_row
  }

  results_df <- .perform_test(object, formula_char, lm_fun, 
                              all_features, assay.type = from)

  # Set a good column order
  variables <- gsub("_P$", "", 
                    colnames(results_df)[grep("P$", colnames(results_df))])
  statistics <- c("Estimate", "LCI95", "UCI95", 
                  "Std_Error", "t_value", "P", "P_FDR")
  col_order <- expand.grid(statistics, variables, stringsAsFactors = FALSE) %>%
    tidyr::unite("Column", "Var2", "Var1")
  col_order <- c("Feature_ID", col_order$Column, c("R2", "Adj_R2"))

  log_text("Linear regression performed.")

  results_df[col_order]
}

#' Linear models ANOVA table
#'
#' Fits a linear model separately for each feature and compute an ANOVA table.
#' Returns all relevant statistics.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param formula_char character, the formula to be used in the linear model 
#' (see Details)
#' @param all_features should all features be included in FDR correction?
#' @param lm_args list of arguments to lm, list names should be parameter names
#' @param anova_args list of arguments to anova, list names should be parameter 
#' names
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A data frame with one row per feature, with all the
#' relevant statistics of the linear model as columns.
#'
#' @details The linear model is fit on combined_data(object). Thus, column names
#' in pheno data can be specified. To make the formulas flexible, the word 
#' "Feature" must be used to signal the role of the features in the formula. 
#' "Feature" will be replaced by the actual Feature IDs during model fitting, 
#' see the example.
#'
#' @examples
#' data(example_set)
#' # A simple example without QC samples
#' # Features predicted by Group and Time
#' lm_anova_results <- perform_lm_anova(drop_qcs(example_set), 
#'   formula_char = "Feature ~ Group + Time")
#'
#' @seealso \code{\link[stats]{lm}}
#'
#' @export
perform_lm_anova <- function(object, formula_char, all_features = FALSE,
                             lm_args = NULL, anova_args = NULL,
                             assay.type = NULL) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from)
  log_text("Starting ANOVA tests")

  anova_fun <- function(feature, formula, data) {
    # Try to fit the linear model
    fit <- NULL
    tryCatch(
      {
        fit <- do.call(stats::lm, 
                       c(list(formula = formula, data = data),
                       lm_args))
        anova_res <- do.call(stats::anova, c(list(object = fit), anova_args))
      },
      error = function(e) message(feature, ": ", e$message))
    if (is.null(anova_res) || sum(!is.na(data[, feature])) < 2) {
      result_row <- NULL
    } else {
      effect_names <- c(names(fit$contrasts), "Residuals")
      anova_names <- c("Df", "Sum_Sq", "Mean_Sq", "F_value", "P")
      result_row <- data.frame(Feature_ID = feature)
      for (i in seq_along(effect_names)) {
        for (j in seq_along(names(anova_res))) {
          name <- paste0(effect_names[i], "_", anova_names[j])
          result_row[name] <- anova_res[[j]][i]
        }
      }
      result_row
    }
    result_row
  }

  results_df <- .perform_test(object, formula_char, anova_fun, 
                              all_features, assay.type = from)

  log_text("ANOVA tests performed")

  results_df
}


#' Logistic regression
#'
#' Fits a logistic regression model separately for each feature. Returns all 
#' relevant statistics.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param formula_char character, the formula to be used in the linear model 
#' (see Details)
#' @param all_features should all features be included in FDR correction?
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional parameters passed to \code{\link[stats]{glm}}
#'
#' @return A data frame with one row per feature, with all the
#' relevant statistics of the linear model as columns.
#'
#' @details The logistic regression model is fit on combined_data(object). 
#' Thus, column names in pheno data can be specified. To make the formulas 
#' flexible, the word "Feature" must be used to signal the role of the features 
#' in the formula. "Feature" will be replaced by the actual Feature IDs during 
#' model fitting, see the example.
#'
#' @examples
#' data(example_set)
#' # A simple example without QC samples
#' # Time predicted by features
#' logistic_results <- perform_logistic(drop_qcs(example_set),
#'   formula_char = "Time ~ Feature + Group"
#' )
#'
#' @seealso \code{\link[stats]{glm}}
#'
#' @export
perform_logistic <- function(object, formula_char, all_features = FALSE, 
                             assay.type = NULL, ...) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from)
  log_text("Starting logistic regression")

  logistic_fun <- function(feature, formula, data) {
    # Try to fit the linear model
    fit <- NULL
    tryCatch(
      {
        fit <- stats::glm(formula, data = data, family = stats::binomial(), ...)
      },
      error = function(e) message(feature, ": ", e$message))
    if (is.null(fit) || sum(!is.na(data[, feature])) < 2) {
      result_row <- NULL
    } else {
      # Gather coefficients and CIs to one data frame row
      coefs <- summary(fit)$coefficients
      confints <- withCallingHandlers(
        expr = stats::confint(fit, level = 0.95, trace = FALSE),
        message = function(m) tryInvokeRestart("muffleMessage"))
      coefs <- data.frame(Variable = rownames(coefs), coefs, 
                          stringsAsFactors = FALSE)
      confints <- data.frame(Variable = rownames(confints), confints,
                             stringsAsFactors = FALSE)

      result_row <- dplyr::left_join(coefs, confints, by = "Variable") %>%
        dplyr::rename("Std_Error" = "Std..Error", "z_value" = "z.value",
                      "P" = "Pr...z..", "LCI95" = "X2.5..", 
                      "UCI95" = "X97.5..") %>%
        tidyr::gather("Metric", "Value", -"Variable") %>%
        tidyr::unite("Column", "Variable", "Metric", sep = "_") %>%
        tidyr::spread("Column", "Value")
      result_row$Feature_ID <- feature
    }
    result_row
  }

  results_df <- .perform_test(object, formula_char, logistic_fun,
                              all_features, assay.type = from)

  # Set a good column order
  variables <- gsub("_P$", "", 
                    colnames(results_df)[grep("P$", colnames(results_df))])
  statistics <- c("Estimate", "LCI95", "UCI95", "Std_Error", 
                  "z_value", "P", "P_FDR")
  col_order <- expand.grid(statistics, variables, stringsAsFactors = FALSE) %>%
    tidyr::unite("Column", "Var2", "Var1")
  col_order <- c("Feature_ID", col_order$Column)
  results_df <- results_df[col_order]
  # Add odds ratios
  estimate_cols <- colnames(results_df)[grepl("_Estimate$",
                                        colnames(results_df))]
  for (estimate_col in estimate_cols) {
    estimate_values <- results_df[, estimate_col]
    results_df <- tibble::add_column(.data = results_df,
                                     OR = exp(estimate_values),
                                     .after = estimate_col)
    estimate_idx <- which(colnames(results_df) == estimate_col)
    colnames(results_df)[estimate_idx + 1] <- gsub("Estimate", "OR",
                                                   estimate_col)
  }

  log_text("Logistic regression performed.")

  results_df
}

.help_lmer <- function(feature, formula, data, ci_method, test_random, ...) {
  # Try to fit the linear model
  fit <- NULL
  # If fitting causes an error, a NULL row is returned
  result_row <- NULL
  tryCatch(
    {
      fit <- lmerTest::lmer(formula, data = data, ...)
    },
    error = function(e) message(feature, ": ", e$message))
  if (!is.null(fit)) {
    # Extract model coefficients
    coefs <- summary(fit)$coefficients
    coefs <- data.frame(Variable = rownames(coefs), coefs, 
                        stringsAsFactors = FALSE)
    # Try to compute confidence intervals
    # If the computation fails, all CIs are NA
    confints <- data.frame(Variable = rownames(coefs), 
                           "X2.5.." = NA, "X97.5.." = NA)
    tryCatch(
      {
        confints <- stats::confint(fit, nsim = 1000, method = ci_method, 
                            oldNames = FALSE)
        confints <- data.frame(Variable = rownames(confints), confints,
                               stringsAsFactors = FALSE)
      },
      error = function(e) message(feature, ": ", e$message)
    )

    # Gather coefficients and CIs to one data frame row
    result_row <- dplyr::left_join(coefs, confints, by = "Variable") %>%
      dplyr::rename("Std_Error" = "Std..Error", "t_value" = "t.value",
                    "P" = "Pr...t..", "LCI95" = "X2.5..", 
                    "UCI95" = "X97.5..") %>%
      tidyr::gather("Metric", "Value", -"Variable") %>%
      tidyr::unite("Column", "Variable", "Metric", sep = "_") %>%
      tidyr::spread("Column", "Value")
    # Add R2 statistics
    result_row$Marginal_R2 <- NA
    result_row$Conditional_R2 <- NA
    tryCatch(
      {
        r2s <- MuMIn::r.squaredGLMM(fit)
        result_row$Marginal_R2 <- r2s[1]
        result_row$Conditional_R2 <- r2s[2]
      },
      error = function(e) message(feature, ": ", e$message)
    )
    # Add Feature ID
    result_row$Feature_ID <- feature
  }

  #Add optional test results for the random effects
  if (test_random) {
    tryCatch(
      {
        r_tests <- as.data.frame(lmerTest::ranova(fit))[-1, c(4, 6)]
        r_tests$Variable <- rownames(r_tests) %>%
          gsub("[(]1 [|] ", "", .) %>%
          gsub("[)]", "", .)
        # Get confidence intervals for the SD of the random effects
        confints$Variable <- confints$Variable %>%
          gsub("sd_[(]Intercept[)][|]", "", .)
        # Get standard deviations of the random effects
        r_variances <- as.data.frame(summary(fit)$varcor)[c("grp", "sdcor")]
        # Join all the information together
        r_result_row <- dplyr::inner_join(r_variances, confints, 
                                          by = c("grp" = "Variable")) %>%
          dplyr::left_join(r_tests, by = c("grp" = "Variable")) %>%
          dplyr::rename(SD = "sdcor", "LCI95" = "X2.5..", "UCI95" = "X97.5..",
                        "P" = "Pr(>Chisq)") %>%
          tidyr::gather("Metric", "Value", "-grp") %>%
          tidyr::unite("Column", "grp", "Metric", sep = "_") %>%
          tidyr::spread("Column", "Value")
        result_row <- cbind(result_row, r_result_row)
      },
      error = function(e) message(feature, ": ", e$message)
    )
  }

  result_row
}

#' Linear mixed models
#'
#' Fits a linear mixed model separately for each feature. Returns all relevant
#' statistics.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param formula_char character, the formula to be used in the linear model 
#' (see Details)
#' @param all_features should all features be included in FDR correction?
#' @param ci_method The method for calculating the confidence intervals as in 
#' \code{\link{confint}}
#' @param test_random logical, whether tests for the significance of the random 
#' effects should be performed
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... additional parameters passed to \code{\link[lmerTest]{lmer}}
#'
#' @return A data frame with one row per feature, with all the
#' relevant statistics of the linear mixed model as columns.
#'
#' @details The model is fit on combined_data(object). Thus, column names
#' in pheno data can be specified. To make the formulas flexible, the word 
#' "Feature" must be used to signal the role of the features in the formula. 
#' "Feature" will be replaced by the actual Feature IDs during model fitting, 
#' see the example. With bootstrap ("boot") confidence intervals, the results 
#' are reproducible if RNGseed is set for the BiocParallel backend. 
#' 
#' @examples
#' data(example_set)
#' # A simple example without QC samples
#' # Features predicted by Group and Time as fixed effects with Subject ID as a 
#' # random effect
#' lmer_results <- perform_lmer(drop_qcs(example_set),
#'   formula_char = "Feature ~ Group + Time + (1 | Subject_ID)",
#'   ci_method = "Wald"
#' )
#' @seealso \code{\link[lmerTest]{lmer}} for model specification 
#'
#' @export
perform_lmer <- function(object, formula_char, all_features = FALSE,
                         ci_method = c("Wald", "profile", "boot"),
                         test_random = FALSE, assay.type = NULL, ...) {
  log_text("Starting fitting linear mixed models.")

  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("Package \'lmerTest\' needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  if (!requireNamespace("MuMIn", quietly = TRUE)) {
    stop("Package \'MuMIn\' needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  .add_citation(paste0("lmerTest package was used for statistical tests in", 
                       " linear mixed models:"), 
                citation("lmerTest"))
  .add_citation(paste0("MuMIn package was used to assess R2 values in",
                       " linear mixed models:"), 
                citation("MuMIn"))

  # Check that ci_method is one of the accepted choices
  ci_method <- match.arg(ci_method)
  
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from)
  results_df <- .perform_test(object, formula_char, .help_lmer, all_features,
                              packages = "lmerTest", ci_method = ci_method,
                              test_random = test_random, assay.type = from)

  # Set a good column order
  fixed_effects <- 
    gsub("_Estimate$", "", 
         colnames(results_df)[grep("Estimate$", colnames(results_df))])
  statistics <- c("Estimate", "LCI95", "UCI95", "Std_Error",
                  "t_value", "P", "P_FDR")
  col_order <- expand.grid(statistics, fixed_effects, 
                           stringsAsFactors = FALSE) %>%
    tidyr::unite("Column", "Var2", "Var1")
  col_order <- c("Feature_ID", col_order$Column, 
                 c("Marginal_R2", "Conditional_R2"))

  if (test_random) {
    random_effects <- 
      gsub("_SD$", "", colnames(results_df)[grep("SD$", colnames(results_df))])
    statistics <- c("SD", "LCI95", "UCI95", "LRT", "P", "P_FDR")
    random_effect_order <- expand.grid(statistics, random_effects,
                                       stringsAsFactors = FALSE) %>%
      tidyr::unite("Column", "Var2", "Var1")
    col_order <- c(col_order, random_effect_order$Column)
  }

  log_text("Linear mixed models fit.")

  results_df[col_order]
}

#' Test homoscedasticity
#'
#' Performs Bartlett's, Levene's and Fligner-Killeen tests for equality of 
#' variances.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param formula_char character, the formula to be used in the linear model 
#' (see Details)
#' @param all_features should all features be included in FDR correction?
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @details The model is fit on combined_data(object). Thus, column names
#' in pheno data can be specified. To make the formulas flexible, the word 
#' "Feature" must be used to signal the role of the features in the formula. 
#' "Feature" will be replaced by the actual Feature IDs during model fitting. 
#' For example, if testing for equality of variances in study groups, use 
#' "Feature ~ Group".
#'
#' @return A data frame with the results.
#'
#' @examples
#' data(example_set)
#' perform_homoscedasticity_tests(example_set, formula_char = "Feature ~ Group")
#'
#' @seealso \code{\link{bartlett.test}}, \code{\link[car]{leveneTest}}, 
#' \code{\link{fligner.test}}
#'
#' @export
perform_homoscedasticity_tests <- function(object, formula_char, 
                                           all_features = FALSE, 
                                           assay.type = NULL) {
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package \'car\' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  .add_citation("car package was used for Levene's test of homoscedasticity:",
                citation("car"))
  log_text("Starting homoscedasticity tests.")
  
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from)

  homosced_fun <- function(feature, formula, data) {
    result_row <- NULL
    tryCatch(
      {
        bartlett <- stats::bartlett.test(formula = formula, data = data)
        levene <- car::leveneTest(y = formula, data = data)
        fligner <- stats::fligner.test(formula = formula, data = data)

        result_row <- data.frame(Feature_ID = feature, 
                                 Bartlett_P = bartlett$p.value,
                                 Levene_P = levene$`Pr(>F)`[1],
                                 Fligner_P = fligner$p.value,
                                 stringsAsFactors = FALSE)
      },
      error = function(e) message(feature, ": ", e$message))
    result_row
  }

  results_df <- .perform_test(object, formula_char, homosced_fun, all_features,
                              assay.type = from)

  log_text("Homoscedasticity tests performed.")

  results_df
}

#' Kruskal-Wallis rank-sum test
#'
#' Performs Kruskal-Wallis rank-sum test for equality.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param formula_char character, the formula to be used in the linear model 
#' (see Details)
#' @param all_features should all features be included in FDR correction?
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @details The model is fit on combined_data(object). Thus, column names
#' in pheno data can be specified. To make the formulas flexible, the word 
#' "Feature" must be used to signal the role of the features in the formula. 
#' "Feature" will be replaced by the actual Feature IDs during model fitting. 
#' For example, if testing for equality of means in study groups, use 
#' "Feature ~ Group".
#'
#' @return A data frame with the results.
#'
#' @seealso \code{\link{kruskal.test}}
#'
#' @examples
#' data(example_set)
#' perform_kruskal_wallis(example_set, formula_char = "Feature ~ Group")
#'
#' @export
perform_kruskal_wallis <- function(object, formula_char, all_features = FALSE,
                                   assay.type = NULL) {
  log_text("Starting Kruskal-Wallis tests.")
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from)

  kruskal_fun <- function(feature, formula, data) {
    result_row <- NULL
    tryCatch(
      {
        kruskal <- stats::kruskal.test(formula = formula, data = data)
        result_row <- data.frame(Feature_ID = feature,
                                 Kruskal_P = kruskal$p.value,
                                 stringsAsFactors = FALSE)
      },
      error = function(e) message(feature, ": ", e$message))

    result_row
  }

  results_df <- .perform_test(object, formula_char, kruskal_fun, all_features,
                              assay.type = from)

  log_text("Kruskal-Wallis tests performed.")

  results_df
}


#' Welch's ANOVA and classic ANOVA
#'
#' Performs ANOVA with Welch's correction as default, to deal with 
#' heterogeneity of variances.
#' Can also perform classic ANOVA with assumption of equal variances.
#' Uses base R function \code{oneway.test}.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param formula_char character, the formula to be used in the linear model 
#' (see Details).
#' @param all_features should all features be included in FDR correction?
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters to \code{\link{oneway.test}}
#'
#' @details The model is fit on combined_data(object). Thus, column names
#' in pheno data can be specified. To make the formulas flexible, the word 
#' "Feature" must be used to signal the role of the features in the formula. 
#' "Feature" will be replaced by the actual Feature IDs during model fitting. 
#' For example, if testing for equality of means in study groups, use 
#' "Feature ~ Group".
#'
#' @return A data frame with the results.
#'
#' @seealso \code{\link{oneway.test}}
#'
#' @examples
#' data(example_set)
#' perform_oneway_anova(example_set, formula_char = "Feature ~ Group")
#'
#' @export
perform_oneway_anova <- function(object, formula_char, all_features = FALSE,
                                 assay.type = NULL, ...) {
  log_text("Starting ANOVA tests.")
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from)

  anova_fun <- function(feature, formula, data) {
    result_row <- NULL
    tryCatch(
      {
        anova_res <- stats::oneway.test(formula = formula, data = data, ...)
        result_row <- data.frame(Feature_ID = feature,
                                 ANOVA_P = anova_res$p.value,
                                 stringsAsFactors = FALSE)
      },
      error = function(e) message(feature, ": ", e$message))

    result_row
  }

  results_df <- .perform_test(object, formula_char, anova_fun, all_features,
                              assay.type = from)

  log_text("ANOVA performed.")

  results_df
}

#' Pairwise and paired t-tests
#'
#' Performs pairwise and paired t-tests. The R default is Welch's t-test 
#' (unequal variances), use var.equal = TRUE for Student's t-test. Use 
#' \code{is_paired} for paired t-tests.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param formula_char character, the formula to be used in the linear model 
#' (see Details)
#' @param is_paired logical, use paired t-test
#' @param id character, name of the subject identification column for paired 
#' version
#' @param all_features should all features be included in FDR correction?
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters passed to \code{\link{t.test}}
#'
#' @details P-values of each comparison are corrected separately from each 
#' other.
#'
#' @return A data frame with the results.
#'
#' @examples
#' data(example_set)
#' # Including QCs as a study group for example
#' t_test_results <- perform_t_test(example_set, 
#'   formula_char = "Feature ~ Group")
#' # Using paired mode (pairs with QC are skipped as there are no common IDs in 
#' # 'example_set')
#' t_test_results <- perform_t_test(example_set,
#'   formula_char = "Feature ~ Time", is_paired = TRUE, id = "Subject_ID")
#' # Only two groups
#' t_test_results <- perform_t_test(drop_qcs(example_set),
#'   formula_char = "Feature ~ Group")
#' 
#' @seealso \code{\link[stats]{t.test}}
#'
#' @export
perform_t_test <- function(object, formula_char, is_paired = FALSE, id = NULL, 
                           all_features = FALSE, assay.type = NULL, ...) {
  message("The functionality of this function has changed.", 
          " It now encompasses pairwise and paired t-tests.")
  message("Remember that t.test returns difference between group means",
          " in different order than lm.\n",
          "This function mimics this behavior, so the effect size is",
          " mean of first level minus mean of second level.")
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from)
  assays(object) <- assays(object)[from]
  
  group <- unlist(strsplit(formula_char, " ~ "))[2]

  if (!is.factor(colData(object)[, group])) {
    stop("Group column should be a factor")
  }

  results_df <- .setup_simple_test(object, group = group, 
                                   id = id, test = "t_test", 
                                   all_features = all_features, 
                                   is_paired = is_paired, ...)
  
  if (length(rownames(object)) != nrow(results_df)) {
    warning("Results don't contain all features.")
  }
  rownames(results_df) <- results_df$Feature_ID
  results_df
}

.adjust_p_values <- function(x, flags) {
  p_cols <- colnames(x)[grepl("_P$", colnames(x))]
  for (p_col in p_cols) {
    p_values <- x[, p_col, drop = TRUE]
    p_values[!is.na(flags)] <- NA
    x <- tibble::add_column(.data = x,
                            FDR = stats::p.adjust(p_values, method = "BH"),
                            .after = p_col)
    p_idx <- which(colnames(x) == p_col)
    colnames(x)[p_idx + 1] <- paste0(p_col, "_FDR")
  }
  x
}

# Calculate pairwise and/or paired statistics
.calc_simple_test <- function(feature, fname, subset1, subset2, 
                              pair, test, is_paired, ...) {
  result_row <- NULL
  tryCatch({
    if (test == "t_test") {
      res <- stats::t.test(feature[subset1], feature[subset2],
                           paired = is_paired, ...)
    } else if (test == "Wilcox") {
      res <- stats::wilcox.test(feature[subset1], feature[subset2],
                                paired = is_paired, conf.int = TRUE, ...)
    }
    # Change name of test for reporting if non-paired non-parametric
    if (test == "Wilcox" & is_paired == FALSE) {
      test <- "Mann_Whitney"
    }
    # Get a single-row data.frame with results
    conf_level <- attr(res$conf.int, "conf.level") * 100
    result_row <- data.frame(Feature_ID = fname, Statistic = res$statistic,
                             Estimate = res$estimate[1], LCI = res$conf.int[1],
                             UCI = res$conf.int[2], P = res$p.value,
                             stringsAsFactors = FALSE)
    ci_idx <- grepl("CI", colnames(result_row))
    colnames(result_row)[ci_idx] <- paste0(colnames(result_row)[ci_idx],
                                           conf_level)
    colnames(result_row)[-1] <- paste0(pair[1], "_vs_", pair[2], "_", 
                                   test, "_", colnames(result_row)[-1])
    result_row
  }, error = function(e) message(fname, ": ", e$message))
} 

# Sets up subsets, calls .calc_simple_test and processes results for simple
# tests
.help_simple_test <- function(object, group, id, test,
                              all_features = FALSE, is_paired, ...) {
  results_df <- NULL
  data <- combined_data(object)
  features <- rownames(object)
  groups <- data[, group]
  pair <- levels(groups)
  if (!is(groups, "factor")) groups <- as.factor(groups)
  # Find samples of complete pairs by group or samples by group
  if (is_paired){
    if (is.null(id)) {
      stop("Subject ID column is missing, please provide it")
    }
    # Get group indices
    subset1 <- which(groups == pair[1])
    subset2 <- which(groups == pair[2])
    # Get complete pairs 
    common_ids <- as.numeric(intersect(data[subset1, id], data[subset2, id]))
    # Keep only complete pairs, order by id
    subset1 <- subset1[data[subset1, id] %in%            
                       common_ids][order(common_ids)]         
    subset2 <- subset2[data[subset2, id] %in% 
                       common_ids][order(common_ids)]         
    log_text(paste0("Starting paired tests for ", 
                    paste0(pair, collapse = " & ")))
    log_text(paste("Found", length(common_ids), "complete pairs."))
    # If there are no complete pairs, return an empty data.frame
    if (length(common_ids) == 0) {
      warning(paste0("Skipped ", paste0(pair, collapse = " & "), 
                     ": no common IDs."))
      return(data.frame("Feature_ID" = features))
    }
  } else {
    # Get group indices
    subset1 <- which(groups == pair[1])
    subset2 <- which(groups == pair[2])
    log_text(paste0("Starting tests for ", paste0(pair, collapse = " & ")))
  }
  
  # Calculate paired or unpaired test for the subsets
  result_rows <- BiocParallel::bpmapply(
    .calc_simple_test, data[, features], features, 
     MoreArgs = list(subset1, subset2, pair, test, is_paired, ...),
     SIMPLIFY = FALSE)
                                          
  results_df <- dplyr::bind_rows(result_rows)
  # Check that results actually contain results
  if (nrow(results_df) == 0) {
     stop("All the tests failed", call. = FALSE)
  }
  rownames(results_df) <- results_df$Feature_ID
  
  # Rows full of NA for features where the test failed
  results_df <- .fill_results(results_df, features)
  # FDR correction
  if (all_features) {
    flags <- rep(NA_character_, nrow(results_df))
  } else {
    flags <- flag(object)
  }
  results_df <- .adjust_p_values(results_df, flags)

  if (is_paired) {
    log_text("Paired tests performed.")
  } else {
    log_text("Tests performed")
  }
  results_df
}

# Calls .help_simple_test in a loop with an object set up for each pairwise
# comparison, also in the case of two groups
.setup_simple_test <- function(object, group, ...) {
  df <- NULL
  groups <- levels(colData(object)[, group])
  combinations <- utils::combn(groups, 2)
  for (i in seq_len(ncol(combinations))) {
    group1 <- as.character(combinations[1, i])
    group2 <- as.character(combinations[2, i])
    # Subset the pair of groups
    object_tmp <- object[, colData(object)[, group] %in% c(group1, group2)]
    colData(object_tmp) <- droplevels(colData(object_tmp))

    res <- .help_simple_test(object_tmp, group = group, ...)
    ifelse(is.null(df), df <- res, df <- dplyr::left_join(df, res))
  }
  df
}

#' Pairwise and paired non-parametric tests
#'
#' Performs pairwise and paired non-parametric tests. The default is Mann-
#' Whitney U test, use \code{is_paired} for Wilcoxon signed rank tests.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or \code{\link{MetaboSet}} object
#' @param formula_char character, the formula to be used in the tests
#' @param is_paired logical, use paired test
#' @param id character, name of the subject identification column for paired 
#' version
#' @param all_features should all features be included in FDR correction?
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters passed to test \code{\link{wilcox.test}}
#'
#' @details P-values of each comparison are corrected separately from each 
#' other. The model is fit on combined_data(object). Thus, column names
#' in pheno data can be specified. To make the formulas flexible, the word 
#' "Feature" must be used to signal the role of the features in the formula. 
#' "Feature" will be replaced by the actual features during model fitting. 
#' For example, if testing for equality of means in study groups, use 
#' "Feature ~ Group".
#'
#' @return A data frame with the results.
#'
#' @examples
#' data(example_set)
#' # Including QCs as a study group for example for pairwise tests
#' mann_whitney_results <- perform_non_parametric(example_set, 
#'   formula_char = "Feature ~ Group")
#' # Using paired mode (pairs with QC are skipped as there are no common IDs in 
#' # 'example_set')
#' wilcoxon_signed_results <- perform_non_parametric(example_set,
#'   formula_char = "Feature ~ Time",
#'   is_paired = TRUE,
#'   id = "Subject_ID")
#' # Only two groups
#' mann_whitney_results <- perform_non_parametric(drop_qcs(example_set), 
#'   formula_char = "Feature ~ Group")
#'
#' @seealso \code{\link[stats]{wilcox.test}}
#'
#' @export
perform_non_parametric <- function(object, formula_char, is_paired = FALSE, 
                                   id = NULL, all_features = FALSE, 
                                   assay.type = NULL, ...) {
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, assay.type = from)
  assays(object) <- assays(object)[from]
  group <- unlist(strsplit(formula_char, " ~ "))[2]

  if (!is.factor(colData(object)[, group])) {
    stop("Grouping column should be a factor")
  }

  results_df <- .setup_simple_test(object, group = group, 
                                   id = id, test = "Wilcox", 
                                   all_features = all_features, 
                                   is_paired = is_paired, ...)
  
  if (length(rownames(object)) != nrow(results_df)) {
    warning("Results don't contain all features.")
  }
  rownames(results_df) <- results_df$Feature_ID
  results_df
}
