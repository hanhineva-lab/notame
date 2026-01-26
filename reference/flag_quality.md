# Flag low-quality features

Flags low-quality features using the quality metrics defined in
(Broadhurst 2018). The metrics are described in more detain in Details.
A condition for keeping the features is given as a character, which is
passed to [`filter`](https://dplyr.tidyverse.org/reference/filter.html).

## Usage

``` r
flag_quality
flag_quality(object, assay.type = NULL, condition =
  "(RSD_r < 0.2 & D_ratio_r < 0.4) | 
  (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)")
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- assay.type:

  character, assay to be used in case of multiple assays

- condition:

  character, condition for keeping the features, see Details

## Value

a SummarizedExperiment object with the features flagged.

## Details

The quality metrics measure two things: internal spread of the QCs, and
spread of the QCs compared to the spread of the biological samples.
Internal spread is measured with relative standard deviation (RSD), also
known as coefficient of variation (CV). \$\$RSD = sd(QC) / mean(QC) \$\$
Where \\sd(QC)\\ is the standard deviation of the QC samples and
'\\mean(QC)\\ is the sample mean of the signal in the QC samples. RSD
can also be replaced by a non-parametric, robust version based on the
median and median absolute deviation (MAD): \$\$RSD_r = 1.4826 \*
MAD(QC) / median(QC)\$\$ The spread of the QC samples compared to the
biological samples is measured using a metric called D-ratio:
\$\$D_ratio = sd(QC) / sd(biological)\$\$ Or, as before, a
non-parametric, robust alternative: \$\$D_ratio_r = MAD(QC) /
MAD(biolofical) \$\$ The default condition keeps features that pass
either of the two following conditions: \$\$RSD_r \< 0.2 \\ D_ratio_r \<
0.4\$\$ \$\$RSD \< 0.1 \\ RSD_r \< 0.1 \\ D_ratio \< 0.1\$\$

## References

Broadhurst, David et al. Guidelines and considerations for the use of
system suitability and quality control samples in mass spectrometry
assays applied in untargeted clinical metabolomic studies. Metabolomics
: Official journal of the Metabolomic Society vol. 14,6 (2018): 72.
doi:10.1007/s11306-018-1367-3

## Examples

``` r
data(toy_notame_set)
ex_set <- flag_quality(toy_notame_set)
#> INFO [2026-01-26 11:40:27] 
#> 92% of features flagged for low quality
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
#> HILIC_neg_259_9623a4_4322 Low_quality  0.127075  0.156465  0.363022  0.475271
#> HILIC_neg_108_1065a2_6121 Low_quality  0.468943  0.375347  1.202952  1.223426
#> HILIC_neg_158_23a1_4128   Low_quality  0.568911  0.188623  1.019238  0.509172
#> HILIC_neg_251_0056a0_6161 Low_quality  0.366231  0.383584  0.875327  1.089626
#> HILIC_neg_401_52a4_211    Low_quality  0.153478  0.222952  0.419890  0.710290
#> ...                               ...       ...       ...       ...       ...
#> RP_pos_226_7218a3_1371    Low_quality  0.419258  0.298690  0.986842  1.069176
#> RP_pos_280_3779a4_0424             NA  0.126026  0.173915  0.253066  0.370787
#> RP_pos_447_3232a2_2393    Low_quality  0.729257  0.948402  1.106946  1.207068
#> RP_pos_453_282a5_128      Low_quality  0.663424  1.133718  1.465057  2.228991
#> RP_pos_153_9636a7_9986    Low_quality  0.570564  0.183891  1.075597  0.439143
# Custom condition
ex_set <- flag_quality(toy_notame_set,
  condition = "RSD_r < 0.3 & D_ratio_r < 0.6")
#> INFO [2026-01-26 11:40:28] 
#> 70% of features flagged for low quality
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
#> HILIC_neg_108_1065a2_6121 Low_quality  0.468943  0.375347  1.202952  1.223426
#> HILIC_neg_158_23a1_4128            NA  0.568911  0.188623  1.019238  0.509172
#> HILIC_neg_251_0056a0_6161 Low_quality  0.366231  0.383584  0.875327  1.089626
#> HILIC_neg_401_52a4_211    Low_quality  0.153478  0.222952  0.419890  0.710290
#> ...                               ...       ...       ...       ...       ...
#> RP_pos_226_7218a3_1371    Low_quality  0.419258  0.298690  0.986842  1.069176
#> RP_pos_280_3779a4_0424             NA  0.126026  0.173915  0.253066  0.370787
#> RP_pos_447_3232a2_2393    Low_quality  0.729257  0.948402  1.106946  1.207068
#> RP_pos_453_282a5_128      Low_quality  0.663424  1.133718  1.465057  2.228991
#> RP_pos_153_9636a7_9986             NA  0.570564  0.183891  1.075597  0.439143
```
