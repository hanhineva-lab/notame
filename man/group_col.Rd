% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constructor.R
\name{group_col}
\alias{group_col}
\alias{group_col<-}
\title{Get and set name of the special column for group labels}
\usage{
group_col(object)

group_col(object) <- value
}
\arguments{
\item{object}{a \code{\link{MetaboSet}} object}

\item{value}{string, name of column to be designated for holding group labels}
}
\value{
Character, the name of the grouping variable.

For the endomorphism, an object with the grouping variable set to 
the specified variable.
}
\description{
Get and set name of the special column for group labels
}
\examples{
data(example_set)
ex_set <- as(example_set, "MetaboSet")
# Get name of grouping variable
group_col(ex_set)

data(example_set)
ex_set <- as(example_set, "MetaboSet")
# Set grouping variable
group_col(ex_set) <- "Group"

}
