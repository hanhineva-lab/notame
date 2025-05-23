% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{impute_rf}
\alias{impute_rf}
\title{Impute missing values using random forest}
\usage{
impute_rf(object, all_features = FALSE, assay.type = NULL, name = NULL, ...)
}
\arguments{
\item{object}{a \code{
\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
or \code{\link{MetaboSet}} object}

\item{all_features}{logical, should all features be used? If FALSE (the 
default), flagged features are removed before imputation.}

\item{assay.type}{character, assay to be used in case of multiple assays}

\item{name}{character, name of the resultant assay in case of multiple assays}

\item{...}{passed to \code{\link[missForest]{missForest}}}
}
\value{
An object as the one supplied, with missing 
values imputed.
}
\description{
Impute the missing values in the peak table of the object using a
random forest. The estimated error in the imputation is logged.
It is recommended to set the seed number for reproducibility
(it is called random forest for a reason).
This a wrapper around \code{\link[missForest]{missForest}}.
Use parallelize = "variables" to run in parallel for faster testing.
NOTE: running in parallel prevents user from setting a seed number.
}
\examples{
data(example_set)
missing <- mark_nas(example_set, 0)
set.seed(38)
imputed <- impute_rf(missing)

}
\seealso{
\code{\link[missForest]{missForest}} for detail about the algorithm 
and the parameters
}
