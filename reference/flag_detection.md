# Flag features with low detection rate

Flags features with too high amount of missing values. There are two
detection rate limits, both defined as the minimum proportion of samples
that need to have a value (not NA) for the feature to be kept.
'`qc_limit` is the detection rate limit for QC samples, '`group_limit`
is the detection rate limit for the actual study groups. If the group
limit is passed for AT LEAST ONE GROUP, then the feature is kept.
Features with low detection rate in QCs are flagged as
"Low_qc_detection", while low detection rate in the study groups is
flagged as "Low_group_detection". The detection rates for all the groups
are recorded in feature data.

## Usage

``` r
flag_detection(
  object,
  qc_limit = 0.7,
  group_limit = 0.5,
  group = NULL,
  assay.type = NULL
)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- qc_limit:

  the detection rate limit for QC samples

- group_limit:

  the detection rate limit for study groups

- group:

  the columns name in sample information to use as the grouping variable

- assay.type:

  character, assay to be used in case of multiple assays

## Value

A SummarizedExperiment object with the features flagged.

## Examples

``` r
data(toy_notame_set)
ex_set <- mark_nas(toy_notame_set, value = 0)
ex_set <- flag_detection(ex_set, group = "Group")
#> INFO [2026-01-26 11:40:27] 
#> 1% of features flagged for low detection rate
rowData(ex_set)
#> DataFrame with 80 rows and 11 columns
#>                                       Feature_ID       Split Alignment
#>                                      <character> <character> <integer>
#> HILIC_neg_259_9623a4_4322 HILIC_neg_259_9623a4..   HILIC_neg         1
#> HILIC_neg_108_1065a2_6121 HILIC_neg_108_1065a2..   HILIC_neg         2
#> HILIC_neg_158_23a1_4128   HILIC_neg_158_23a1_4..   HILIC_neg         3
#> HILIC_neg_251_0056a0_6161 HILIC_neg_251_0056a0..   HILIC_neg         4
#> HILIC_neg_401_52a4_211    HILIC_neg_401_52a4_211   HILIC_neg         5
#> ...                                          ...         ...       ...
#> RP_pos_226_7218a3_1371    RP_pos_226_7218a3_1371      RP_pos        16
#> RP_pos_280_3779a4_0424    RP_pos_280_3779a4_0424      RP_pos        17
#> RP_pos_447_3232a2_2393    RP_pos_447_3232a2_2393      RP_pos        18
#> RP_pos_453_282a5_128        RP_pos_453_282a5_128      RP_pos        19
#> RP_pos_153_9636a7_9986    RP_pos_153_9636a7_9986      RP_pos        20
#>                           Average_Mz Average_Rt_min      Column    Ion_mode
#>                            <numeric>      <numeric> <character> <character>
#> HILIC_neg_259_9623a4_4322    259.962       4.432241       HILIC         neg
#> HILIC_neg_108_1065a2_6121    108.107       2.612073       HILIC         neg
#> HILIC_neg_158_23a1_4128      158.230       1.412822       HILIC         neg
#> HILIC_neg_251_0056a0_6161    251.006       0.616145       HILIC         neg
#> HILIC_neg_401_52a4_211       401.520       4.210970       HILIC         neg
#> ...                              ...            ...         ...         ...
#> RP_pos_226_7218a3_1371       226.722        3.13714          RP         pos
#> RP_pos_280_3779a4_0424       280.378        4.04238          RP         pos
#> RP_pos_447_3232a2_2393       447.323        2.23926          RP         pos
#> RP_pos_453_282a5_128         453.282        5.12797          RP         pos
#> RP_pos_153_9636a7_9986       153.964        7.99864          RP         pos
#>                                  Flag Detection_rate_Group_A
#>                           <character>              <numeric>
#> HILIC_neg_259_9623a4_4322          NA                   0.95
#> HILIC_neg_108_1065a2_6121          NA                   0.90
#> HILIC_neg_158_23a1_4128            NA                   0.95
#> HILIC_neg_251_0056a0_6161          NA                   0.90
#> HILIC_neg_401_52a4_211             NA                   0.95
#> ...                               ...                    ...
#> RP_pos_226_7218a3_1371             NA                   0.85
#> RP_pos_280_3779a4_0424             NA                   0.85
#> RP_pos_447_3232a2_2393             NA                   0.95
#> RP_pos_453_282a5_128               NA                   0.90
#> RP_pos_153_9636a7_9986             NA                   0.80
#>                           Detection_rate_Group_B Detection_rate_QC
#>                                        <numeric>         <numeric>
#> HILIC_neg_259_9623a4_4322                   0.90               1.0
#> HILIC_neg_108_1065a2_6121                   0.95               0.9
#> HILIC_neg_158_23a1_4128                     0.90               0.8
#> HILIC_neg_251_0056a0_6161                   0.95               1.0
#> HILIC_neg_401_52a4_211                      0.95               1.0
#> ...                                          ...               ...
#> RP_pos_226_7218a3_1371                       0.9               0.9
#> RP_pos_280_3779a4_0424                       0.9               1.0
#> RP_pos_447_3232a2_2393                       0.9               0.9
#> RP_pos_453_282a5_128                         1.0               0.8
#> RP_pos_153_9636a7_9986                       0.9               0.8
```
