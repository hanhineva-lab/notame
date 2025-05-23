% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constructor.R
\name{read_from_excel}
\alias{read_from_excel}
\title{Read formatted Excel files}
\usage{
read_from_excel(
  file,
  sheet = 1,
  id_column = NULL,
  corner_row = NULL,
  corner_column = NULL,
  id_prefix = "ID_",
  split_by = NULL,
  name = NULL,
  mz_limits = c(10, 2000),
  rt_limits = c(0, 20),
  skip_checks = FALSE
)
}
\arguments{
\item{file}{path to the Excel file}

\item{sheet}{the sheet number or name}

\item{id_column}{character, column name for unique identification of samples}

\item{corner_row}{integer, the bottom row of sample information,
usually contains data file names and feature info column names.
If set to NULL, will be detected automatically.}

\item{corner_column}{integer or character, the corresponding column number 
or the column name (letter) in Excel.
If set to NULL, will be detected automatically.}

\item{id_prefix}{character, prefix for autogenerated sample IDs, see Details}

\item{split_by}{character vector, in the case where all the modes are in the 
same Excel file, the column names of feature data used to separate the modes 
(usually Mode and Column)}

\item{name}{in the case where the Excel file only contains one mode, the 
name of the mode, such as "Hilic_neg"}

\item{mz_limits}{numeric vector of two, all m/z values should be in between 
these}

\item{rt_limits}{numeric vector of two, all retention time values should be 
in between these}

\item{skip_checks}{logical: skip checking and fixing of data integrity. Not 
recommended, but sometimes useful when you just want to read the data in as 
is and fix errors later. The data integrity checks are important for 
functioning of notame.}
}
\value{
A list of three data frames:
\itemize{
  \item exprs: the actual abundances, size features x samples
  \item pheno_data: sample information, size sample info variables x samples
  \item feature_data: information about the features, size features x 
feature info columns
}
}
\description{
Reads data from an Excel file of the following format:
\itemize{
  \item Left side of the sheet contains information about the features, size 
features x feature info columns
  \item Top part contains sample information, size sample info variables x 
samples
  \item The middle contains the actual abundances, size features x samples
}
This function separates the three parts from the file, and returns them in a 
list.
}
\details{
If skip_checks = FALSE, \code{\link{read_from_excel}} attempts to modify the 
data as per \code{\link{fix_object}} and checks the data. If skip_checks 
= TRUE, parameters for \code{fix_object} are ignored.
}
\examples{
data <- read_from_excel(
  file = system.file("extdata", "example_set.xlsx", 
  package = "notame"), sheet = 1, corner_row = 11, corner_column = "H",
  split_by = c("Column", "Ion_mode"))

modes <- construct_metabosets(exprs = data$exprs, 
  pheno_data = data$pheno_data, feature_data = data$feature_data,
  group_col = "Group")

}
