# Finish a log

Logs the current date and time and session info, and switches logging
off.

## Usage

``` r
finish_log()
```

## Value

None, the function is invoked for its side effect.

## See also

[`init_log`](https://hanhineva-lab.github.io/notame/reference/init_log.md),
[`log_text`](https://hanhineva-lab.github.io/notame/reference/log_text.md)

## Examples

``` r
finish_log()
#> INFO [2026-01-26 11:40:26] Finished analysis. Mon Jan 26 11:40:26 2026
#> Session info:
#> INFO [2026-01-26 11:40:26] R Under development (unstable) (2026-01-25 r89330)
#> INFO [2026-01-26 11:40:26] Platform: x86_64-pc-linux-gnu
#> INFO [2026-01-26 11:40:26] Running under: Ubuntu 24.04.3 LTS
#> INFO [2026-01-26 11:40:26] 
#> INFO [2026-01-26 11:40:26] Matrix products: default
#> INFO [2026-01-26 11:40:26] BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> INFO [2026-01-26 11:40:26] LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> INFO [2026-01-26 11:40:26] 
#> INFO [2026-01-26 11:40:26] locale:
#> INFO [2026-01-26 11:40:26]  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#> INFO [2026-01-26 11:40:26]  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=C              
#> INFO [2026-01-26 11:40:26]  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#> INFO [2026-01-26 11:40:26]  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#> INFO [2026-01-26 11:40:26]  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> INFO [2026-01-26 11:40:26] [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> INFO [2026-01-26 11:40:26] 
#> INFO [2026-01-26 11:40:26] time zone: UTC
#> INFO [2026-01-26 11:40:26] tzcode source: system (glibc)
#> INFO [2026-01-26 11:40:26] 
#> INFO [2026-01-26 11:40:26] attached base packages:
#> INFO [2026-01-26 11:40:26] [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> INFO [2026-01-26 11:40:26] [8] base     
#> INFO [2026-01-26 11:40:26] 
#> INFO [2026-01-26 11:40:26] other attached packages:
#> INFO [2026-01-26 11:40:26]  [1] notame_1.1.3                SummarizedExperiment_1.41.0
#> INFO [2026-01-26 11:40:26]  [3] Biobase_2.71.0              GenomicRanges_1.63.1       
#> INFO [2026-01-26 11:40:26]  [5] Seqinfo_1.1.0               IRanges_2.45.0             
#> INFO [2026-01-26 11:40:26]  [7] S4Vectors_0.49.0            BiocGenerics_0.57.0        
#> INFO [2026-01-26 11:40:26]  [9] generics_0.1.4              MatrixGenerics_1.23.0      
#> INFO [2026-01-26 11:40:26] [11] matrixStats_1.5.0           ggplot2_4.0.1              
#> INFO [2026-01-26 11:40:26] 
#> INFO [2026-01-26 11:40:26] loaded via a namespace (and not attached):
#> INFO [2026-01-26 11:40:26]  [1] gtable_0.3.6         xfun_0.56            bslib_0.10.0        
#> INFO [2026-01-26 11:40:26]  [4] httr2_1.2.2          htmlwidgets_1.6.4    lattice_0.22-7      
#> INFO [2026-01-26 11:40:26]  [7] notameStats_1.1.1    vctrs_0.7.1          tools_4.6.0         
#> INFO [2026-01-26 11:40:26] [10] parallel_4.6.0       curl_7.0.0           fansi_1.0.7         
#> INFO [2026-01-26 11:40:26] [13] tibble_3.3.1         pkgconfig_2.0.3      Matrix_1.7-4        
#> INFO [2026-01-26 11:40:26] [16] RColorBrewer_1.1-3   S7_0.2.1             desc_1.4.3          
#> INFO [2026-01-26 11:40:26] [19] lifecycle_1.0.5      compiler_4.6.0       farver_2.1.2        
#> INFO [2026-01-26 11:40:26] [22] textshaping_1.0.4    fontawesome_0.5.3    codetools_0.2-20    
#> INFO [2026-01-26 11:40:26] [25] htmltools_0.5.9      sass_0.4.10          yaml_2.3.12         
#> INFO [2026-01-26 11:40:26] [28] pillar_1.11.1        pkgdown_2.2.0        jquerylib_0.1.4     
#> INFO [2026-01-26 11:40:26] [31] whisker_0.4.1        BiocParallel_1.45.0  openssl_2.3.4       
#> INFO [2026-01-26 11:40:26] [34] DelayedArray_0.37.0  cachem_1.1.0         abind_1.4-8         
#> INFO [2026-01-26 11:40:26] [37] notameViz_1.1.1      tidyselect_1.2.1     digest_0.6.39       
#> INFO [2026-01-26 11:40:26] [40] dplyr_1.1.4          purrr_1.2.1          fastmap_1.2.0       
#> INFO [2026-01-26 11:40:26] [43] grid_4.6.0           cli_3.6.5            SparseArray_1.11.10 
#> INFO [2026-01-26 11:40:26] [46] magrittr_2.0.4       S4Arrays_1.11.1      withr_3.0.2         
#> INFO [2026-01-26 11:40:26] [49] scales_1.4.0         rappdirs_0.3.4       lambda.r_1.2.4      
#> INFO [2026-01-26 11:40:26] [52] rmarkdown_2.30       XVector_0.51.0       igraph_2.2.1        
#> INFO [2026-01-26 11:40:26] [55] otel_0.2.0           futile.logger_1.4.9  ragg_1.5.0          
#> INFO [2026-01-26 11:40:26] [58] askpass_1.2.1        memoise_2.0.1        evaluate_1.0.5      
#> INFO [2026-01-26 11:40:26] [61] knitr_1.51           viridisLite_0.4.2    futile.options_1.0.1
#> INFO [2026-01-26 11:40:26] [64] rlang_1.1.7          downlit_0.4.5        glue_1.8.0          
#> INFO [2026-01-26 11:40:26] [67] formatR_1.14         BiocManager_1.30.27  xml2_1.5.2          
#> INFO [2026-01-26 11:40:26] [70] rstudioapi_0.18.0    jsonlite_2.0.0       R6_2.6.1            
#> INFO [2026-01-26 11:40:26] [73] systemfonts_1.3.1    fs_1.6.6            
```
