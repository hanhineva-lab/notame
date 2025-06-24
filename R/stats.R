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
#' object
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