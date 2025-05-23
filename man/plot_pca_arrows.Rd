% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scatter.R
\name{plot_pca_arrows}
\alias{plot_pca_arrows}
\title{PCA plot with arrows}
\usage{
plot_pca_arrows(
  object,
  pcs = c(1, 2),
  all_features = FALSE,
  center = TRUE,
  scale = "uv",
  color,
  time,
  subject,
  alpha = 0.6,
  arrow_style = arrow(),
  title = "PCA changes",
  subtitle = NULL,
  color_scale = getOption("notame.color_scale_dis"),
  text_base_size = 14,
  line_width = 0.5,
  assay.type = NULL,
  ...
)
}
\arguments{
\item{object}{a \code{
\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
or \code{\link{MetaboSet}} object}

\item{pcs}{numeric vector of length 2, the principal components to plot}

\item{all_features}{logical, should all features be used? If FALSE (the 
default), flagged features are removed before visualization.}

\item{center}{logical, should the data be centered prior to PCA? (usually 
yes)}

\item{scale}{scaling used, as in \code{\link[pcaMethods]{prep}}. Default is 
"uv" for unit variance}

\item{color}{character, name of the column used for coloring the arrows}

\item{time}{character, name of the column containing timepoints}

\item{subject}{character, name of the column containing subject identifiers}

\item{alpha}{numeric, value for the alpha parameter of the arrows 
(transparency)}

\item{arrow_style}{a description of arrow heads, the size and angle can be 
modified, see \code{?arrow}}

\item{title, subtitle}{the titles of the plot}

\item{color_scale}{the color scale as returned by a ggplot function}

\item{text_base_size}{the base size of the text}

\item{line_width}{the width of the arrows}

\item{assay.type}{character, assay to be used in case of multiple assays}

\item{...}{additional arguments passed to \code{\link[pcaMethods]{pca}}}
}
\value{
A ggplot object.
}
\description{
Plots changes in PCA space according to time. All the observations of a 
single subject are connected by an arrow ending at the last observation.
}
\examples{
data(example_set)
plot_pca_arrows(drop_qcs(example_set), color = "Group", time = "Time",
  subject = "Subject_ID")
# If the sample size is large, plot groups separately
plot_pca_arrows(drop_qcs(example_set), color = "Group", 
                time = "Time", subject = "Subject_ID") +
  facet_wrap(~Group)

}
\seealso{
\code{\link[pcaMethods]{pca}}
}
