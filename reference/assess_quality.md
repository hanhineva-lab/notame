# Assess quality information of features

Assess features using the quality metrics defined in (Broadhurst 2018).
The quality metrics are described in Details section of
[`flag_quality`](https://hanhineva-lab.github.io/notame/reference/flag_quality.md)

## Usage

``` r
assess_quality(object, assay.type = NULL)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- assay.type:

  character, assay to be used in case of multiple assays

## Value

A SummarizedExperiment object with quality metrics in feature data.

## Examples

``` r
data(toy_notame_set)
ex_set <- assess_quality(toy_notame_set)
rowData(ex_set)
#> DataFrame with 80 rows and 12 columns
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
#>                                  Flag       RSD     RSD_r   D_ratio D_ratio_r
#>                           <character> <numeric> <numeric> <numeric> <numeric>
#> HILIC_neg_259_9623a4_4322          NA  0.127075  0.156465  0.363022  0.475271
#> HILIC_neg_108_1065a2_6121          NA  0.468943  0.375347  1.202952  1.223426
#> HILIC_neg_158_23a1_4128            NA  0.568911  0.188623  1.019238  0.509172
#> HILIC_neg_251_0056a0_6161          NA  0.366231  0.383584  0.875327  1.089626
#> HILIC_neg_401_52a4_211             NA  0.153478  0.222952  0.419890  0.710290
#> ...                               ...       ...       ...       ...       ...
#> RP_pos_226_7218a3_1371             NA  0.419258  0.298690  0.986842  1.069176
#> RP_pos_280_3779a4_0424             NA  0.126026  0.173915  0.253066  0.370787
#> RP_pos_447_3232a2_2393             NA  0.729257  0.948402  1.106946  1.207068
#> RP_pos_453_282a5_128               NA  0.663424  1.133718  1.465057  2.228991
#> RP_pos_153_9636a7_9986             NA  0.570564  0.183891  1.075597  0.439143
```
