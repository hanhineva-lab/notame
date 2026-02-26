# Flag contaminants based on blanks

Flags contaminant features by comparing either median, mean or max of
blanks and biological samples. Biological samples are defined as samples
that are not marked as blanks and are not QCs.

## Usage

``` r
flag_contaminants(
  object,
  blank_col,
  blank_label,
  blank_type = c("mean", "median", "max"),
  sample_type = c("max", "median", "mean"),
  flag_thresh = 5,
  flag_label = "Contaminant",
  assay.type = NULL
)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- blank_col:

  character, the column name in colData with blank labels

- blank_label:

  character, the label for blank samples in blank_col

- blank_type:

  character, one of "mean", "median", or "max"

- sample_type:

  character, one of "max", "median", or "mean"

- flag_thresh:

  numeric, the scaled ratio threshold for flagging contaminants. If
  blank_type(blanks) \* flag_thresh \> sample_type(biological samples),
  the feature gets flagged.

- flag_label:

  character, the label used when flagging contaminants. Can be changed
  if sample processing contaminants and carryover contaminants are
  flagged separately.

- assay.type:

  character, assay to be used in case of multiple assays

## Value

A SummarizedExperiment object with contaminant features flagged.

## Details

If the calculation(biological samples) \< the calculation(blanks) times
a set ratio, the feature is flagged as contaminant. Default calculations
are "max" for biological samples and "mean" for blanks.

## Examples

``` r
data(toy_notame_set)
# Make a blank sample which has one (first) feature exceeding the threshold
## Abundance matrix
max_sample <- max(assay(toy_notame_set)[1, toy_notame_set$QC != "QC"])
assay <- matrix(c(max_sample * 0.20 + 1, rep(0, 79)), ncol = 1, nrow = 80, 
                  dimnames = list(NULL, "Demo_51"))
assay <- cbind(assay(toy_notame_set), assay)
## Sample metadata
pheno_data <- colData(toy_notame_set)[1, ]
rownames(pheno_data) <- "Demo_51"
pheno_data$Sample_ID <- "Demo_51"
pheno_data$Injection_order <- 51
pheno_data[c("Subject_ID", "Group", "QC", "Time")] <- "Blank"
pheno_data <- rbind(colData(toy_notame_set), pheno_data)
## Feature metadata
feature_data <- rowData(toy_notame_set)

# Construct SummarizedExperiment object with blank sample
ex_set <- SummarizedExperiment(assays = assay, 
                               colData = pheno_data,
                               rowData = feature_data)
# Flag contaminant(s)
contaminants_flagged <- flag_contaminants(ex_set, blank_col = "QC", 
                                          blank_label = "Blank")
#> INFO [2026-02-26 10:48:33] 
#> 1% of features flagged as contaminants
```
