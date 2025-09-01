#' Toy data set
#'
#' Contains imaginary data used in testing the package functions. The dataset 
#' has 50 samples and 80 features. This dataset includes multiple observations 
#' from same subjects, sampled at two timepoints in separate batches and 
#' divided to two groups. The analytical modes are also available as separate 
#' SummarizedExperiment objects.
#' Note that across batches, the features don't have different feature ID's, 
#' m/z and retention time as would be the case with real-world data. In 
#' essence, the example data reflects that features were aligned perfectly 
#' between batches. 
"toy_notame_set"

#' @rdname toy_notame_set
"hilic_neg_sample"

#' @rdname toy_notame_set
"hilic_pos_sample"

#' @rdname toy_notame_set
"rp_neg_sample"

#' @rdname toy_notame_set
"rp_pos_sample"
