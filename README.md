# NOTE This is a development version. There is active development with breaking changes until the package has been approved in Bioconductor. Yet, everything in the main branch is to our knowledge working as it should. The original notame package this development is based on can be installed using `devtools::install_github("antonvsdata/notame@v0.3.1")`

# notame - Workflow for non-targeted LC-MS metabolic profiling 

This package can be used to analyze preprocessed LC-MS data in non-targeted metabolomics. Notame was developed at the [research group of nutritional metabolomics at University of Eastern Finland](https://www3.uef.fi/en/web/kttravi/metabolomics2) and [Afekta Technologies](https://afekta.com/), a spinoff metabolomics company. We use notame as a way to bundle together all the preprocessing methods we use for our non-targeted LC-MS metabolomics data, so it mainly consists of methods found in other packages, and a bunch of visualizations we have found useful.

For more detailed information on how we run our LC-MS experiments and where this package fits in our workflow, you can find the paper here: ["notame": Workflow for Non-Targeted LC-MS Metabolic Profiling](https://www.mdpi.com/2218-1989/10/4/135). A huge thank you for everyone involved in the paper!

Currently, the package is developed by Afekta Technologies. The package API is still quite experimental, and breaking changes are possible, although new features are added quite slowly since the package fulfills its current tasks well.

### What does notame do?

Before we go into the list of features, it is good for you to know hot the workflow in our lab works. The first step is to take raw data files created by the LC-MS instrument and create a peak table using a peak picking software (we use [MS-DIAL](http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/)). After peak picking with the dedicated software, we use R for data preprocessing, quality control, statistical analysis and visualization. We then use the obtained results in identification of the actual metabolites. During the years, we ended up with various scripts that were hard to handle and update, so we decided to make this package to keep things under control. 

Here is a list of the current main functionalities of the package:

- Reading data from Excel spreadsheets created with MS-DIAL
- Data is stored in a custom object that holds all the information about the features and samples along with the feature abundance matrix. This allows for a simple interface for all of the functions in the package, as there is no need to juggle with different matrices/data frames.
- Drift correction: correcting for systematic drift in the intensity of molecular features using cubic spline correction (see [Kirwan & Broadhurst et al.](https://doi.org/10.1007/s00216-013-6856-7))
- Identifying and flagging (or removing) low-quality molecular features using quality metrics defined by [Broadhurst et al.](https://doi.org/10.1007/s11306-018-1367-3)
- Imputing missing values, multiple strategies available. Random forest imputation recommended, see [Kokla et al.](https://doi.org/10.1186/s12859-019-3110-0)
- Batch effect correction: correcting for systematic variation between batches. Multiple strategies available.
- A novel method for clustering similar molecular features
- A bunch of statistical analyses, both feature-wise tests and multivariate models
- A rather nice set of visualizations for use in quality control, explorative analysis and interpretation of results from statistical tests


## Installation and getting started

PACKAGE requires R version 3.5.0 or greater.

notame functions depend on a ton of other R packages. The packages you need to install depend on what you're using notame for: some packages are only needed for specific visualizations, others for batch effect correction, and others for common preprocessing tasks. This is why it's recommended to only install the packages that are really needed to make notame work. To do this, run:

```
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("hanhineva-lab/notame")
```

After installing the package, you can install rest of the packages you need on the fly OR use a handy function called ```install_dependencies```, which lets you install packages for core preprocessing, batch correction, specific visualizations etc.

### Installing development version

To install the current development version between releases, install the package from the dev branch:

```
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("hanhineva-lab/notame", ref = "dev")
```

## Documentation

For instructions on how to use the package, run `browseVignettes("notame")`. The Introduction vignette should get you started pretty well!


## Credits and license

The first version of notame was written by Anton Klåvus for his master's thesis in Bioinformatics at Aalto university (published under former name Anton Mattsson), while working for University of Eastern Finland and Afekta Technologies. Notame is inspired by analysis scripts written by Jussi Paananen, Oskari Timonen and Anton Klåvus at University of Eastern Finland. The algorithm for clustering molecular features originating from the same compound is based on MATLAB code written by David Broadhurst, Professor of Data Science & Biostatistics in the School of Science, and director of the Centre for Integrative Metabolomics & Computational Biology at the Edith Covan University.

If you find any bugs or other things to fix, please submit an issue on GitHub! All contributions to the package are always welcome!

notame is published under an MIT license (tl;dr: it's really permissive!)


