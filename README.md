# NOTE This is a development version. There is active development with breaking changes until the package has been approved in Bioconductor. The original notame package this development is based on can be installed using `devtools::install_github('antonvsdata/notame@v0.3.1')`

# notame - Workflow for non-targeted LC-MS metabolic profiling 

This package can be used to analyze preprocessed LC-MS data in non-targeted metabolomics. Notame was developed at the [research group of nutritional metabolomics at University of Eastern Finland](https://hanhinevalab.com/home) and [Afekta Technologies](https://afekta.com/), a spinoff metabolomics company. We use notame as a way to bundle together all the preprocessing methods we use for our non-targeted LC-MS metabolomics data, so it mainly consists of methods found in other packages, and a bunch of visualizations we have found useful.

For more detailed information on how we run our LC-MS experiments and where this package fits in our workflow, you can find the paper here: ["notame": Workflow for Non-Targeted LC-MS Metabolic Profiling](https://www.mdpi.com/2218-1989/10/4/135). A huge thank you for everyone involved in the paper!

Currently, the package is developed by Afekta Technologies and was reworked for Bioconductor to promote interoperability with functionality that can complement the notame workflow. This resulted in some breaking changes (see NEWS for user-facing changes). See the vignettes and documentation for more information. 

### What does notame do?
Before we go into the list of features, it is good for you to know how the workflow in our lab works. The first step is to take raw data files created by the LC-MS instrument and create a peak table using a peak picking software (we use [MS-DIAL](https://systemsomicslab.github.io/compms/msdial/main.html)). After peak picking with the dedicated software, we use R for data preprocessing, quality control, statistical analysis and visualization. We then use the obtained results in identification of the actual metabolites. During the years, we ended up with various scripts that were hard to handle and update, so we decided to make this package to keep things under control. 

Here is a list of the current main functionalities of the package:

- Reading data from Excel spreadsheets created with MS-DIAL
- Data is stored in a SummarizedExperiment object that holds all the information about the features and samples along with the feature abundance matrix. This allows for a simple interface for all of the functions in the package, as there is no need to juggle with different matrices/data frames.
- Drift correction: correcting for systematic drift in the intensity of molecular features using cubic spline correction (see [Kirwan & Broadhurst et al.](https://doi.org/10.1007/s00216-013-6856-7))
- Identifying and flagging (or removing) low-quality molecular features using quality metrics defined by [Broadhurst et al.](https://doi.org/10.1007/s11306-018-1367-3)
- Imputing missing values, multiple strategies available. Random forest imputation recommended, see [Kokla et al.](https://doi.org/10.1186/s12859-019-3110-0)
- Batch effect correction: correcting for systematic variation between batches
- A novel method for clustering similar molecular features
- A bunch of statistical analyses, both feature-wise tests and multivariate models
- A rather nice set of visualizations for use in quality control, explorative analysis and interpretation of results from statistical tests


## Installation and getting started

### Bioc-release

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("notame")
```

### Bioc-devel

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("notame")
```

## Credits and license

The first version of notame was written by Anton Klåvus for his master's thesis in Bioinformatics at Aalto university (published under former name Anton Mattsson), while working for University of Eastern Finland and Afekta Technologies. Notame is inspired by analysis scripts written by Jussi Paananen, Oskari Timonen and Anton Klåvus at University of Eastern Finland. The algorithm for clustering molecular features originating from the same compound is based on MATLAB code written by David Broadhurst. Development was picked up by Atte Lihtamo, Retu Haikonen, Vilhelm Suksi and Leo Lahti, leading to the Bioconductor release. 

If you find any bugs or other things to fix, please submit an issue on GitHub! All contributions to the package are always welcome!

notame is published under an MIT license (tl;dr: it's really permissive!)


