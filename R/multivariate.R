# ----------- Random Forest ----------------

#' Fit Random Forest
#'
#' Fits a random forest, where given response column in pheno data is predicted 
#' using the features. Can be used both for classification and regression. For 
#' more information, see the documentation of \code{randomForest::randomForest}.
#' After fitting the random forest, use rf_importance as a shortcut for getting 
#' the feature importance in random forest prediction.
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param y character, column name of phenoData giving the dependent variable 
#' of the model
#' @param all_features logical, should all features be included in the model? 
#' if FALSE, flagged features are left out
#' @param covariates character, column names of pheno data
#' to use as covariates in the model, in addition to molecular features
#' @param importance Should importance of features be assessed?
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters passed to \code{randomForest::randomForest}
#'
#' @return An object of class randomForest.
#'
#' @seealso \code{\link[randomForest]{randomForest}}, 
#' \code{\link{importance_rf}}
#'
#' @examples
#' rf <- fit_rf(example_set, y = "Group")
#' rf
#' importance_rf(rf)
#'
#' @export
fit_rf <- function(object, y, all_features = FALSE, 
                   covariates = NULL, importance = TRUE, 
                   assay.type = NULL, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package \"randomForest\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("randomForest package was used to fit random forest models:",
                citation("randomForest"))
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(y, covariates),
                         assay.type = from)
  object <- drop_flagged(object, all_features = all_features)

  x <- combined_data(object, assay.type = from)
  x <- x[, c(rownames(object), covariates)]
  rf <- randomForest::randomForest(x = x, y = colData(object)[, y], 
                                   importance = importance, ...)

  rf
}


#' Feature importance in random forest
#'
#' Extracts feature importance in random forest in a nice format.
#'
#' @param rf An object of class randomForest
#'
#' @return A data frame of feature importance.
#'
#' @seealso \code{\link[randomForest]{randomForest}}, \code{\link{fit_rf}}
#'
#' @examples
#' rf <- fit_rf(example_set, y = "Group")
#' rf
#' importance_rf(rf)
#'
#' @export
importance_rf <- function(rf) {
  # Extract metrics and feature ID
  df <- data.frame(Feature_ID = rownames(rf$importance),
                   as.data.frame(rf$importance),
                   stringsAsFactors = FALSE, check.names = FALSE)
  df
}


# ------------------ mixOmics PLS ---------------------------

#' A helper function for extracting predictor matrix with covariates
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param covariates character, column names of pheno data to use as covariates 
#' in the model, in addition to molecular features
#' @return A data frame with predictors, including covariates.
#' @noRd
.get_x <- function(object, covariates, assay.type = NULL) {
  # Convert covariates to numeric
  if (any(!vapply(colData(object)[, covariates], .looks_numeric, logical(1)))) {
    stop("All covariates should be convertable to numeric.")
  }
  colData(object)[covariates] <- lapply(colData(object)[covariates], as.numeric)

  # Extract X
  x <- combined_data(object, assay.type)[, c(rownames(object), covariates)]
  x
}

#' Plot points in PLS space
#'
#' A helper function for \code{mixomics_pls} and \code{mixomics_spls}.
#'
#' @param model a PLS or sPLS model
#' @param Y the Y matrix
#' @param y the name of the y variable
#' @param title plot title
#' @noRd
.plot_pls <- function(model, Y, y, title) {
  # Extract scores and add y variable
  scores <- data.frame(model$variates$X[, seq_len(2)])
  colnames(scores) <- c("X1", "X2")
  scores[, y[1]] <- Y[, 1]
  # Explained variance as percentage
  var_exp <- 100 * model$prop_expl_var$X[seq_len(2)] %>% round(digits = 3)
  p <- ggplot(scores, aes(x = .data[["X1"]], y = .data[["X2"]], 
                          color = .data[[y]])) +
    geom_point() +
    getOption("notame.color_scale_con") +
    theme_minimal() +
    labs(x = paste("X1:", var_exp[1], "%"),
         y = paste("X2:", var_exp[2], "%"),
         title = title)
  grDevices::dev.new()
  plot(p)
}


#' PLS
#'
#' Simple wrappers for fitting a PLS model using pls function of the mixOmics 
#' package. The object can then be passed to many of the mixOmics functions for 
#' prediction, performance evaluation etc. Also plot a scores plot of the first 
#' two components.
#' \itemize{
#' \item{\code{mixomics_pls} A simple PLS model with set number of components 
#' and all features}
#' \item{\code{mixomics_pls_optimize} Test different numbers of components}
#' \item{\code{mixomics_spls_optimize} sPLS model: Test different numbers of 
#' components and features}
#' }
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param y character vector, column names of the grouping variable to predict
#' @param ncomp number of X components
#' @param folds the number of folds to use in k-fold cross validation
#' @param nrepeat the number of times to repeat the cross validation. Lower 
#' this for faster testing.
#' @param plot_scores logical, if TRUE, a scatter plot with the first two PLS-
#' components as x and y-axis will be drawn, colored by the y-variable. 
#' Only really makes sense if y is a single variable
#' @param all_features logical, should all features be included in the model? 
#' if FALSE, flagged features are left out
#' @param covariates character, column names of pheno datato use as covariates 
#' in the model, in addition to molecular features
#' @param n_features the number of features to try for each component
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... any parameters passed to \code{mixOmics::pls} or 
#' \code{mixOmics::spls}
#'
#' @return An object of class "mixo_pls" or "mixo_spls".
#'
#' @examples
#' pls_res <- mixomics_pls(example_set, y = "Injection_order", ncomp = 3)
#' pls_opt <- mixomics_pls_optimize(example_set, 
#'   y = "Injection_order", ncomp = 3)
#' pls_res <- mixomics_spls_optimize(example_set,
#'   y = "Injection_order", ncomp = 3,
#'   n_features <- c(1:10, 12, 15, 20)
#' )
#' @name pls
#' @seealso \code{\link[mixOmics]{pls}}, \code{\link[mixOmics]{perf}},
#' \code{\link[mixOmics]{spls}}, \code{\link[mixOmics]{tune.spls}}
NULL

#' @rdname pls
#' @export
mixomics_pls <- function(object, y, ncomp, plot_scores = TRUE, 
                         all_features = FALSE, covariates = NULL, 
                         assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))

  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                         assay.type = from)

  predictors <- .get_x(object, covariates, from)
  outcome <- colData(object)[y]

  log_text("Fitting PLS")
  pls_model <- mixOmics::pls(predictors, outcome, ncomp = ncomp, ...)

  if (plot_scores && ncomp > 1) {
    .plot_pls(pls_model, outcome, y, 
              title = "PLS: first 2 components and the outcome variable")
  }
  pls_model
}

#' @rdname pls
#'
#' @export
mixomics_pls_optimize <- function(object, y, ncomp, folds = 5, nrepeat = 50,
                                  plot_scores = TRUE, all_features = FALSE,
                                  covariates = NULL, assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                         assay.type = from)
  pls_res <- mixomics_pls(object = object, y = y, ncomp = ncomp, 
                          plot_scores = FALSE, all_features = all_features,
                          covariates = covariates, assay.type = from, ...)

  log_text("Evaluating PLS performance")
  perf_pls <- mixOmics::perf(pls_res, validation = "Mfold", 
                             folds = folds, nrepeat = nrepeat)

  # Plot Mean Square Error
  p1 <- ggplot(data.frame(ncomp = seq_len(ncomp),
                          MSEP = as.vector(perf_pls$measure$MSEP$summary$mean)),
               aes(x = ncomp, y = .data$MSEP)) +
    geom_line() +
    labs(color = NULL, title = "Mean Square Error") +
    theme_bw() +
    scale_x_continuous(breaks = seq_len(ncomp)) +
    theme(panel.grid.minor.x = element_blank())

  # Plot R2 and Q2
  plot_data <- data.frame(R2 = as.vector(perf_pls$measure$R2$summary$mean),
                          Q2 = as.vector(perf_pls$measure$Q2$summary$mean),
                          ncomp = seq_len(ncomp)) %>%
    tidyr::gather(key = "key", value = "value", -ncomp)

  p2 <- ggplot(plot_data, aes(x = ncomp, y = .data$value, color = .data$key)) +
    geom_line() +
    labs(color = NULL, title = "R2 and Q2") +
    theme_bw() +
    getOption("notame.color_scale_dis") +
    scale_x_continuous(breaks = seq_len(ncomp)) +
    theme(panel.grid.minor.x = element_blank())
  
  plot(cowplot::plot_grid(p1, p2, nrow = 1))

  # Find the optimal number of components
  ncomp_opt <- which(perf_pls$measure$MSEP$summary$mean ==
    min(perf_pls$measure$MSEP$summary$mean))[1]
  log_text(paste0("Choosing a PLS model with ", ncomp_opt, 
                  " component(s) based on the minimal MSE\n",
                  "Take a look at the plot and make sure this", 
                  " is the correct number of components"
  ))
  
  mixomics_pls(object = object, y = y, ncomp = ncomp_opt, 
               plot_scores = plot_scores, covariates = covariates, ...)
}

#' @rdname pls
#'
#' @export
mixomics_spls_optimize <- function(object, y, ncomp, n_features =
                                   c(seq_len(10), seq(20, 300, 10)), 
                                   folds = 5, nrepeat = 50, plot_scores = TRUE,
                                   all_features = FALSE, covariates = NULL,
                                   assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                         assay.type = from)

  predictors <- .get_x(object, covariates, from)
  outcome <- colData(object)[y]

  # Test different number of components and features with cross validation
  log_text("Tuning sPLS")
  tuned_spls <- mixOmics::tune.spls(predictors, outcome, ncomp = ncomp,
                                    test.keepX = n_features,
                                    validation = "Mfold", folds = folds,
                                    nrepeat = nrepeat, measure = "MAE")
  # Plot error for each component with different number of features
  plot(plot(tuned_spls) + ggtitle("Performance of sPLS models"))
  # Choose optimal numbers of components and features
  ncomp_opt <- tuned_spls$choice.ncomp$ncomp
  keep_x <- tuned_spls$choice.keepX[seq_len(ncomp_opt)]
  log_text(paste("Final model has", ncomp_opt, 
                 "components with the numbers of features:",
                 paste(keep_x, collapse = ", ")))
  # Fit the final model
  spls_final <- mixOmics::spls(predictors, outcome, ncomp = ncomp_opt, 
                               keepX = keep_x, ...)
  # Scatter plot of points in PLS space
  if (ncomp_opt > 1) {
    .plot_pls(spls_final, outcome, y, title = 
             "Final sPLS model: first 2 components and the outcome variable")
  }

  spls_final
}

.plot_plsda <- function(model, y, title, dist = "max.dist") {
  background <- mixOmics::background.predict(model, comp.predicted = 2, 
                                             dist = dist)
  grDevices::dev.new()
  mixOmics::plotIndiv(model, comp = seq_len(2), group = y, 
                      ind.names = FALSE, title = paste(title), 
                      legend = TRUE, ellipse = TRUE)
  grDevices::dev.new()
  mixOmics::plotIndiv(model, comp = seq_len(2), group = y, ind.names = FALSE,
                      title = paste(title, "prediction areas"), 
                      legend = TRUE, background = background)
}

#' PLS-DA
#'
#' A simple wrapper for fitting a PLS-DA model using plsda function of the 
#' mixOmics package. The object can then be passed to many of the mixOmics 
#' functions for prediction, performance evaluation etc.
#' \itemize{
#' \item{\code{mixomics_plsda} A simple PLS-DA model with set number of 
#' components and all features}
#' \item{\code{mixomics_plsda_optimize} Test different numbers of components,
#' choose the one with minimal balanced error rate}
#' \item{\code{mixomics_splsda_optimize} Test different numbers 
#' of components and features, choose the one with minimal balanced error rate}
#' }
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param y character, column name of the grouping variable to predict
#' @param ncomp the number of X components
#' @param folds the number of folds to use in k-fold cross validation
#' @param nrepeat the number of times to repeat the cross validation. Lower 
#' this for faster testing.
#' @param n_features the number of features to try for each component
#' @param dist the distance metric to use, one of "max.dist", 
#' "mahalanobis.dist", "centroids.dist".
#' use \code{\link{mixomics_plsda_optimize}} to find the best distance metric
#' @param plot_scores logical, if TRUE, a scatter plot with the first two PLS-
#' components as x and y-axis will be drawn, with both prediction surface and 
#' ellipses
#' @param all_features logical, should all features be included in the model? 
#' if FALSE, flagged features are left out
#' @param covariates character, column names of pheno data to use as covariates 
#' in the model, in addition to molecular features
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... any parameters passed to \code{mixOmics::plsda}
#'
#' @return An object of class "mixo_plsda".
#'
#' @examples
#' noqc <- drop_qcs(example_set)
#' plsda_res <- mixomics_plsda(noqc, y = "Group", ncomp = 2)
#' set.seed(38)
#' plsda_opt <- mixomics_plsda_optimize(noqc, y = "Group", ncomp = 3)
#' set.seed(38)
#' splsda_opt <- mixomics_splsda_optimize(noqc,
#'   y = "Group", dist = "max.dist", ncomp = 2,
#'   n_features <- c(1:10, 12, 15, 20)
#' )
#' 
#' @name pls_da
#' @seealso \code{\link[mixOmics]{plsda}}, \code{\link[mixOmics]{perf}},
#' \code{\link[mixOmics]{splsda}}, \code{\link[mixOmics]{tune.splsda}}
NULL

#' @rdname pls_da
#' @export
mixomics_plsda <- function(object, y, ncomp, plot_scores = TRUE, 
                           all_features = FALSE, covariates = NULL,
                           assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.",
         "Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)  
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                         assay.type = from)

  predictors <- .get_x(object, covariates, from)
  outcome <- colData(object)[, y]
  # outcome needs to be a factor, this ensures the levels are right
  if (!is(outcome, "factor")) {
    outcome <- as.factor(outcome)
    warning(y, " is not encoded as a factor, converted to factor with levels: ",
            paste(levels(outcome), collapse = ", "))
  }
  log_text("Fitting PLS-DA")
  plsda_model <- mixOmics::plsda(predictors, outcome, ncomp = ncomp, ...)
  if (plot_scores && ncomp > 1) {
    .plot_plsda(plsda_model, outcome, title = "PLS-DA")
  }

  plsda_model
}


#' @rdname pls_da
#' @export
mixomics_plsda_optimize <- function(object, y, ncomp, folds = 5, nrepeat = 50,
                                    plot_scores = TRUE, all_features = FALSE,
                                    covariates = NULL, assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)  
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                         assay.type = from)

  plsda_res <- mixomics_plsda(object = object, y = y, ncomp = ncomp,
                              plot_scores = FALSE, all_features = all_features,
                              covariates = covariates, 
                              assay.type = from, ...)

  log_text("Evaluating PLS-DA performance")
  perf_plsda <- mixOmics::perf(plsda_res, validation = "Mfold", folds = 5,
                               auc = TRUE, nrepeat = 50)

  plot(perf_plsda, col = mixOmics::color.mixo(seq_len(3)), 
       sd = TRUE, legend.position = "horizontal")
  graphics::title("Performance of PLS-DA models")
  # Find the distance metric with minimum BER
  ber <- perf_plsda$error.rate$BER
  inds <- which(ber == min(ber), arr.ind = TRUE)[1, ]
  dist_met <- colnames(ber)[inds[2]]
  # Find the optimal number of components
  ncomp_opt <- perf_plsda$choice.ncomp["BER", dist_met]
  log_text(paste("Choosing a PLS-DA model with", ncomp_opt, 
                 "components using", dist_met))

  mixomics_plsda(object = object, y = y, ncomp = ncomp_opt, 
                 plot_scores = plot_scores, assay.type = from, ...)
}


#' @rdname pls_da
#' @export
mixomics_splsda_optimize <- function(object, y, ncomp, dist,
                                     n_features = c(seq_len(10), 
                                                    seq(20, 300, 10)),
                                     folds = 5, nrepeat = 50,
                                     plot_scores = TRUE, all_features = FALSE,
                                     covariates = NULL, assay.type = NULL, 
                                     ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object,  pheno_cols = c(y, covariates), 
                         assay.type = from)

  predictors <- .get_x(object, covariates, from)
  outcome <- colData(object)[, y]
  # outcome needs to be a factor, this ensures the levels are right
  if (!is(outcome, "factor")) {
    outcome <- as.factor(outcome)
    warning(y, " is not encoded as a factor, converted to factor with levels: ",
            paste(levels(outcome), collapse = ", "))
  }
  # Test different components and numbers of features with cross validation
  log_text("Tuning sPLS-DA")
  tuned_splsda <- mixOmics::tune.splsda(predictors, outcome, ncomp = ncomp,
                                        validation = "Mfold", folds = folds,
                                        dist = dist, measure = "BER", 
                                        nrepeat = nrepeat,
                                        test.keepX = n_features)
  # Plot error rate of different components as a function of number of features
  plot(plot(tuned_splsda) + ggtitle("Performance of sPLS-DA models"))
  
  # Choose optimal numbers of components and features
  ncomp_opt <- tuned_splsda$choice.ncomp$ncomp
  keep_x <- tuned_splsda$choice.keepX[seq_len(ncomp_opt)]
  log_text(paste("Final model has", ncomp_opt, 
                 "components with the numbers of features:",
                 paste(keep_x, collapse = ", ")))
  # Fit the final model
  splsda_final <- mixOmics::splsda(predictors, outcome, ncomp = ncomp_opt,
                                   keepX = keep_x)
  # Scatterplots with prediction surface and ellipses
  if (plot_scores && ncomp_opt > 1) {
    grDevices::dev.new()
    .plot_plsda(splsda_final, outcome, title = "Final sPLS-DA model", 
               dist = dist)
  }
  splsda_final
}



# ------------------- MUVR --------------------------------

#' Multivariate modelling with minimally biased variable selection (MUVR)
#'
#' A wrapper around the MUVR2 (random forest, PLS(-DA)) and MUVR2_EN (elastic 
#' net) functions from the MUVR2 package. 
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param y character, column name in pheno data of the target variable
#' @param id character, column name in pheno data of the subject ID variable in 
#' case of repeated measurements
#' @param multi_level logical, whether multi-level modeling should be applied, 
#' see Details
#' @param multi_level_var character, column name in pheno data of the variable 
#' for splitting the data in multi-level modeling
#' @param all_features logical, should all features be included in the model? 
#' if FALSE, flagged features are left out
#' @param covariates,static_covariates character, column names of 
#' pheno data to use as covariates in the model, in addition to molecular 
#' features. 
#' \code{static_covariates} are ignored for non-multi-level models.
#' For multi-level models, the change in \code{covariates} is computed, while
#' \code{static_covariates} are taken from the first time point. 
#' @param nRep Number of repetitions of double CV, parameter of MUVR
#' @param nOuter Number of outer CV loop segments, parameter of MUVR
#' @param nInner Number of inner CV loop segments, parameter of MUVR
#' @param varRatio Ratio of variables to include in subsequent inner loop 
#' iteration, parameter of MUVR
#' @param method Multivariate method. Supports 'PLS', 'RF' and 'EN'
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters to MUVR2::MUVR2 or MUVR2::MUVR2_EN 
#' and MUVR2::getVar (when method == "EN")
#'
#' @return A MUVR object. (make this more descriptive)
#'
#' @details This function is now using the MUVR2 package, characterized as an 
#' upgrade extending the original MUVR package by the inclusion of elastic net 
#' regression (EN) and other functionality not covered by this wrapper. Elastic 
#' net regression supports covariate adjustment by suppressing regularization 
#' of specified features from the regularization procedure. Note that this is 
#' different from simply including covariates such as sex. EN also differs from 
#' PLS and RF in that no recursive variable elimination is performed, so an 
#' additional scheme is used to obtain the 'min', 'mid' and 'max' models using
#' MUVR2::getVar in this wrapper. 
#' 
#' Sex would be entered as a static covariate, since the change in sex is zero 
#' for all individuals, so computing the change and using that as a covariate 
#' does not make sense.
#'
#' Note that there are several more plots available in MUVR2 for inspecting the 
#' results, notably MUVR2::plotMV, MUVR2::plotStability and MUVR2::plotVIRank. 
#' Many of these return different plots depending on the model specification.
#'
#' @examples
#' # PLS simple
#' rf_model <- muvr_analysis(drop_qcs(example_set), 
#'   y = "Group", nRep = 2, method = "PLS")
#'
#' # RF with covariate and repeated measures (not longitudinal)
#' ex_set <- drop_qcs(example_set)
#' ex_set$Injection_order %<>% as.numeric()
#' en_model <- muvr_analysis(drop_qcs(ex_set), y = "Group", id = "Subject_ID", 
#'   nRep = 2, method = "RF", covariates = "Injection_order")
#'
#' # RF on multilevel variable comparing levels of y
#' rf_model_ <- muvr_analysis(drop_qcs(example_set), 
#'   y = "Group", multi_level = TRUE, id = "Subject_ID", 
#'   multi_level_var = "Time", method = "RF", nRep = 2)
#'
#' # EN on multilevel variable with covariate and static covariate
#' ex_set <- drop_qcs(example_set)
#' ex_set$Injection_order %<>% as.numeric()
#' ex_set$Group %<>% as.numeric()
#' en_model <- muvr_analysis(ex_set, id = "Subject_ID", 
#'  multi_level = TRUE, multi_level_var = "Time", 
#'  covariates = "Injection_order", static_covariates = "Group", 
#'  method = "EN", nRep = 2)
#'
#' @seealso \code{\link[MUVR2]{MUVR2}} \code{\link[MUVR2]{MUVR2_EN}} 
#' \code{\link[MUVR2]{getVar}} \code{\link[MUVR2]{plotMV}} 
#' \code{\link[MUVR2]{plotStability}} \code{\link[MUVR2]{plotVIRank}}
#'
#' @export
muvr_analysis <- function(object, y = NULL, id = NULL, multi_level = FALSE,
                          multi_level_var = NULL, covariates = NULL,
                          static_covariates = NULL, all_features = FALSE,
                          nRep = 50, nOuter = 6, nInner = nOuter - 1,
                          varRatio = 0.75, method = c("PLS", "RF"),
                          assay.type = NULL, ...) {
  if (!requireNamespace("MUVR2", quietly = TRUE)) {
    stop("Package \"MUVR2\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation(paste("MUVR2 package was used to fit multivariate models",
                      "with variable selection:"),
                citation("MUVR2"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, 
                         pheno_cols = c(y, id, covariates, static_covariates,
                                        multi_level_var, covariates), 
                         assay.type = from)

  # MUVR2 can only use numeric predictors
  classes <- vapply(colData(object)[, c(covariates, static_covariates)], 
                    class, character(1))
  if (length(classes) && any(classes != "numeric")) {
    stop("MUVR2 can only deal with numeric inputs,", 
         " please transform all covariates to numeric", call. = FALSE)
  }

  if (any(!vapply(colData(object)[, covariates], .looks_numeric, logical(1)))) {
    stop("All covariates should be convertable to numeric.")
  }
  colData(object)[covariates] <- lapply(colData(object)[covariates], as.numeric)
  
  # Do do.call with MUVR2::MUVR2_EN if method == "EN", to avoid nesting
  if (method == "EN") {
    func <- MUVR2::MUVR2_EN
  } else {
    func <- MUVR2::MUVR2
  }
  
  # Substitute additional arguments to formal arguments of MUVR2 or MUVR2_EN
  add_args <- substitute(...())
  add_args <- add_args[names(add_args) %in% formalArgs(func)]
  
  # Classic MUVR
  if (!multi_level) {
    if (is.null(y)) {
      stop("y variable needs to be defined unless doing multi-level modeling.")
    }
    predictors <- combined_data(object, from)[, c(rownames(object), covariates)]
    outcome <- colData(object)[, y]

    # Independent samples
    if (is.null(id)) {
      # Make list of arguments for do.call, match additional arguments to 
      # formal arguments of MUVR2 or MUVR2_EN
      args <- list(X = predictors, Y = outcome, nRep = nRep, nOuter = nOuter, 
                   nInner = nInner, varRatio = varRatio, method = method)
      muvr_model <- do.call(func, c(args, add_args))
    } else {
      # Multiple measurements
      ID <- as.numeric(colData(object)[, id])
      args <- list(X = predictors, Y = outcome, ID = ID, nRep = nRep, 
                   nOuter = nOuter, nInner = nInner, varRatio = varRatio, 
                   method = method)

      muvr_model <- do.call(func, c(args, add_args))
    }
  } else { # Multi-level analysis
    if (is.null(id) || is.null(multi_level_var)) {
      stop("id and multi_level_var needed for multi-level modeling.")
    }
    # Check that multi_level_var has only 2 unique values
    ml_var <- colData(object)[, multi_level_var] <-
      as.factor(colData(object)[, multi_level_var])
    if (length(levels(ml_var)) != 2) {
      stop("The multilevel variable should have exactly 2 unique values.")
    } else {
      message("Computing effect matrix according to ", multi_level_var, " : ",
              levels(ml_var)[2], " - ", levels(ml_var)[1])
    }

    # Compute effect matrix with covariates
    cd <- combined_data(object, from)
    cd <- cd[order(cd[, id]), ]
    x1 <- cd[cd[, multi_level_var] == levels(ml_var)[1],
             c(rownames(object), covariates)]
    x2 <- cd[cd[, multi_level_var] == levels(ml_var)[2],
             c(rownames(object), covariates)]
    predictors <- x2 - x1
    # Add static covariates, where we don't want to compute change, such as sex
    predictors[, static_covariates] <- 
      cd[cd[, multi_level_var] == levels(ml_var)[1], static_covariates]
    rownames(predictors) <- unique(cd[, id])

    # Modeling
    if (!is.null(y)) { # Compare change of multi_level_var between levels of y
      outcome <- cd[cd[, multi_level_var] == levels(ml_var)[1], y]
      args <- list(X = predictors, Y = outcome, nRep = nRep, nOuter = nOuter,
                   nInner = nInner, varRatio = varRatio, method = method)
      muvr_model <- do.call(func, c(args, add_args))
    } else { # Compare levels of multi_level_var
      args <- list(X = predictors, ML = TRUE, nRep = nRep, nOuter = nOuter,
                   nInner = nInner, varRatio = varRatio, method = method)
      muvr_model <- do.call(func, c(args, add_args))
    }
  }
  
  if (method == "EN") {
    add_args <- substitute(...())
    add_args <- add_args[names(add_args) %in% formalArgs(MUVR2::getVar)]    
    muvr_model <- do.call(MUVR2::getVar, 
                          c(list(rdCVnetObject = muvr_model), add_args))
  }
  
  # Plot performance
  MUVR2::plotVAL(muvr_model)
  
  muvr_model
}

# PERMANOVA ----

#' PERMANOVA
#'
#' Performs permutational multivariate analysis of variance. Uses package 
#' called PERMANOVA.
#'
#' @param object a SummarizedExperiment or MetaboSet object
#' @param group character, name of the column to compare
#' @param all_features should all features be included?
#' @param transform Transformation to use in \code{IniTransform}. By default 
#' uses "Standardize columns".
#' @param coef Coefficient to calculate continuous distances in 
#' \code{DistContinuous}.
#' By default uses Pythagorean distances.
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters to \code{\link[PERMANOVA]{PERMANOVA}}
#'
#' @return A PERMANOVA object.
#'
#' @examples
#' permanova_res <- perform_permanova(drop_qcs(example_set), group = "Group")
#'
#' @export
perform_permanova <- function(object, group, all_features = FALSE,
                              transform = "Standardize columns",
                              coef = "Pythagorean", assay.type = NULL, ...) {
  if (!requireNamespace("PERMANOVA", quietly = TRUE)) {
    stop("Package \"PERMANOVA\" needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  
  .add_citation(paste("PERMANOVA was used for permutational multivariate", 
                      "analysis of variance:"),
                citation("PERMANOVA"))
                
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_factors = group,
                         assay.type = from)

  log_text("Starting PERMANOVA tests")
  data <- t(assay(object, from))
  data <- PERMANOVA::IniTransform(data, transform = transform)
  initialized <- PERMANOVA::DistContinuous(data, coef = coef)
  res <- PERMANOVA::PERMANOVA(initialized, colData(object)[, group], ...)
  log_text("PERMANOVA performed")

  res
}
