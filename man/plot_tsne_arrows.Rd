% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scatter.R
\name{plot_tsne_arrows}
\alias{plot_tsne_arrows}
\title{t-SNE plot with arrows}
\usage{
plot_tsne_arrows(
  object,
  all_features = FALSE,
  center = TRUE,
  scale = "uv",
  perplexity = 30,
  pca_method = "nipals",
  color,
  time,
  subject,
  alpha = 0.6,
  arrow_style = arrow(),
  title = "t-SNE changes",
  subtitle = paste("Perplexity:", perplexity),
  color_scale = getOption("notame.color_scale_dis"),
  text_base_size = 14,
  line_width = 0.5,
  assay.type = NULL,
  ...
)
}
\arguments{
\item{object}{a SummarizedExperiment or MetaboSet object}

\item{all_features}{logical, should all features be used? If FALSE (the 
default), flagged features are removed before visualization.}

\item{center}{logical, should the data be centered prior to PCA? (usually 
yes)}

\item{scale}{scaling used, as in  \code{\link[pcaMethods]{prep}}. Default is 
"uv" for unit variance}

\item{perplexity}{the perplexity used in t-SNE}

\item{pca_method}{the method used in PCA if there are missing values}

\item{color}{character, name of the column used for coloring the points}

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

\item{...}{additional arguments passed to \code{\link[Rtsne]{Rtsne}}}
}
\value{
A ggplot object. If \code{density} is \code{TRUE}, the plot will 
consist of multiple parts and is harder to modify.
}
\description{
Computes t-SNE into two dimensions and plots changes according to time.
All the observations of a single subject are connected by an arrow ending at 
the last observation. In case there are missing values, PCA is performed 
using the nipals method of \code{\link[pcaMethods]{pca}}, the method can be 
changed to "ppca" if nipals fails.
}
\examples{
data(example_set)
plot_tsne_arrows(drop_qcs(example_set), perplexity = 10, color = "Group", 
  time = "Time", subject = "Subject_ID")
# If the sample size is large, plot groups separately
plot_tsne_arrows(drop_qcs(example_set), perplexity = 10, color = "Group", 
  time = "Time", subject = "Subject_ID") +
    facet_wrap(~Group)

}
\seealso{
\code{\link[Rtsne]{Rtsne}}
}
