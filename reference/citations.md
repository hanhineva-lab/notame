# Show citations

This function lists citations behind the notame functions that have been
called during the session. All notame functions update the list
automatically. The citations are taken from the call to
'`citation("package")`, and complemented with a brief description of
what the package was used for. NOTE: the citations might not point to
the correct paper if the package authors have not supplied correct
citation information for their package. The output is written to the
current log file, if specified.

## Usage

``` r
citations()
```

## Value

None, the function is invoked for its side effect.

## Examples

``` r
citations()
#> INFO [2026-02-26 10:48:28] Preprocessing and analyses were performed using notame package:
#> INFO [2026-02-26 10:48:28] To cite package ‘notame’ in publications use:
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   Klavus, A.; Kokla, M.; Noerman, S.; Koistinen, V.M.; Tuomainen, M.;
#> INFO [2026-02-26 10:48:28]   Zarei, I.; Meuronen, T.; Hakkinen, M.R.; Rummukainen, S.; Farizah
#> INFO [2026-02-26 10:48:28]   Babu, A.; Sallinen, T.; Karkkainen, O.; Paananen, J.; Broadhurst, D.;
#> INFO [2026-02-26 10:48:28]   Brunius, C.; Hanhineva, K. "notame": Workflow for Non-Targeted LC-MS
#> INFO [2026-02-26 10:48:28]   Metabolic Profiling. Metabolites 2020, 10, 135.
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28] A BibTeX entry for LaTeX users is
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   @Article{,
#> INFO [2026-02-26 10:48:28]     title = {"notame": Workflow for Non-Targeted LC-MS Metabolic Profiling},
#> INFO [2026-02-26 10:48:28]     author = {Anton Klavus and Marietta Kokla and Stefania Noerman and Ville M. Koistinen and Marjo Tuomainen and Iman Zarei and Topi Meuronen and Merja R. Hakkinen and Soile Rummukainen and Ambrin {Farizah Babu} and Taisa Sallinen and Olli Karkkainen and Jussi Paananen and David Broadhurst and Carl Brunius and Kati Hanhineva},
#> INFO [2026-02-26 10:48:28]     year = {2020},
#> INFO [2026-02-26 10:48:28]     journal = {Metabolites},
#> INFO [2026-02-26 10:48:28]     volume = {10},
#> INFO [2026-02-26 10:48:28]     number = {135},
#> INFO [2026-02-26 10:48:28]     url = {https://doi.org/10.3390/metabo10040135},
#> INFO [2026-02-26 10:48:28]     doi = {10.3390/metabo10040135},
#> INFO [2026-02-26 10:48:28]   }
#> INFO [2026-02-26 10:48:28] The primary data structure in notame is SummarizedExperiment:
#> INFO [2026-02-26 10:48:28] To cite package ‘SummarizedExperiment’ in publications use:
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   Morgan M, Obenchain V, Hester J, Pagès H (2026).
#> INFO [2026-02-26 10:48:28]   _SummarizedExperiment: A container (S4 class) for matrix-like
#> INFO [2026-02-26 10:48:28]   assays_. doi:10.18129/B9.bioc.SummarizedExperiment
#> INFO [2026-02-26 10:48:28]   <https://doi.org/10.18129/B9.bioc.SummarizedExperiment>. R package
#> INFO [2026-02-26 10:48:28]   version 1.41.1,
#> INFO [2026-02-26 10:48:28]   <https://bioconductor.org/packages/SummarizedExperiment>.
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28] A BibTeX entry for LaTeX users is
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   @Manual{,
#> INFO [2026-02-26 10:48:28]     title = {SummarizedExperiment: A container (S4 class) for matrix-like assays},
#> INFO [2026-02-26 10:48:28]     author = {Martin Morgan and Valerie Obenchain and Jim Hester and Hervé Pagès},
#> INFO [2026-02-26 10:48:28]     year = {2026},
#> INFO [2026-02-26 10:48:28]     note = {R package version 1.41.1},
#> INFO [2026-02-26 10:48:28]     url = {https://bioconductor.org/packages/SummarizedExperiment},
#> INFO [2026-02-26 10:48:28]     doi = {10.18129/B9.bioc.SummarizedExperiment},
#> INFO [2026-02-26 10:48:28]   }
#> INFO [2026-02-26 10:48:28] visualizations in notame are built with ggplot2:
#> INFO [2026-02-26 10:48:28] To cite ggplot2 in publications, please use
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
#> INFO [2026-02-26 10:48:28]   Springer-Verlag New York, 2016.
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28] A BibTeX entry for LaTeX users is
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   @Book{,
#> INFO [2026-02-26 10:48:28]     author = {Hadley Wickham},
#> INFO [2026-02-26 10:48:28]     title = {ggplot2: Elegant Graphics for Data Analysis},
#> INFO [2026-02-26 10:48:28]     publisher = {Springer-Verlag New York},
#> INFO [2026-02-26 10:48:28]     year = {2016},
#> INFO [2026-02-26 10:48:28]     isbn = {978-3-319-24277-4},
#> INFO [2026-02-26 10:48:28]     url = {https://ggplot2.tidyverse.org},
#> INFO [2026-02-26 10:48:28]   }
data(toy_notame_set)
ex_set <- flag_quality(toy_notame_set)
#> INFO [2026-02-26 10:48:28] 
#> 92% of features flagged for low quality
# Broadhurst et al.(2018) added to citations
citations()
#> INFO [2026-02-26 10:48:28] Preprocessing and analyses were performed using notame package:
#> INFO [2026-02-26 10:48:28] To cite package ‘notame’ in publications use:
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   Klavus, A.; Kokla, M.; Noerman, S.; Koistinen, V.M.; Tuomainen, M.;
#> INFO [2026-02-26 10:48:28]   Zarei, I.; Meuronen, T.; Hakkinen, M.R.; Rummukainen, S.; Farizah
#> INFO [2026-02-26 10:48:28]   Babu, A.; Sallinen, T.; Karkkainen, O.; Paananen, J.; Broadhurst, D.;
#> INFO [2026-02-26 10:48:28]   Brunius, C.; Hanhineva, K. "notame": Workflow for Non-Targeted LC-MS
#> INFO [2026-02-26 10:48:28]   Metabolic Profiling. Metabolites 2020, 10, 135.
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28] A BibTeX entry for LaTeX users is
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   @Article{,
#> INFO [2026-02-26 10:48:28]     title = {"notame": Workflow for Non-Targeted LC-MS Metabolic Profiling},
#> INFO [2026-02-26 10:48:28]     author = {Anton Klavus and Marietta Kokla and Stefania Noerman and Ville M. Koistinen and Marjo Tuomainen and Iman Zarei and Topi Meuronen and Merja R. Hakkinen and Soile Rummukainen and Ambrin {Farizah Babu} and Taisa Sallinen and Olli Karkkainen and Jussi Paananen and David Broadhurst and Carl Brunius and Kati Hanhineva},
#> INFO [2026-02-26 10:48:28]     year = {2020},
#> INFO [2026-02-26 10:48:28]     journal = {Metabolites},
#> INFO [2026-02-26 10:48:28]     volume = {10},
#> INFO [2026-02-26 10:48:28]     number = {135},
#> INFO [2026-02-26 10:48:28]     url = {https://doi.org/10.3390/metabo10040135},
#> INFO [2026-02-26 10:48:28]     doi = {10.3390/metabo10040135},
#> INFO [2026-02-26 10:48:28]   }
#> INFO [2026-02-26 10:48:28] The primary data structure in notame is SummarizedExperiment:
#> INFO [2026-02-26 10:48:28] To cite package ‘SummarizedExperiment’ in publications use:
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   Morgan M, Obenchain V, Hester J, Pagès H (2026).
#> INFO [2026-02-26 10:48:28]   _SummarizedExperiment: A container (S4 class) for matrix-like
#> INFO [2026-02-26 10:48:28]   assays_. doi:10.18129/B9.bioc.SummarizedExperiment
#> INFO [2026-02-26 10:48:28]   <https://doi.org/10.18129/B9.bioc.SummarizedExperiment>. R package
#> INFO [2026-02-26 10:48:28]   version 1.41.1,
#> INFO [2026-02-26 10:48:28]   <https://bioconductor.org/packages/SummarizedExperiment>.
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28] A BibTeX entry for LaTeX users is
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   @Manual{,
#> INFO [2026-02-26 10:48:28]     title = {SummarizedExperiment: A container (S4 class) for matrix-like assays},
#> INFO [2026-02-26 10:48:28]     author = {Martin Morgan and Valerie Obenchain and Jim Hester and Hervé Pagès},
#> INFO [2026-02-26 10:48:28]     year = {2026},
#> INFO [2026-02-26 10:48:28]     note = {R package version 1.41.1},
#> INFO [2026-02-26 10:48:28]     url = {https://bioconductor.org/packages/SummarizedExperiment},
#> INFO [2026-02-26 10:48:28]     doi = {10.18129/B9.bioc.SummarizedExperiment},
#> INFO [2026-02-26 10:48:28]   }
#> INFO [2026-02-26 10:48:28] visualizations in notame are built with ggplot2:
#> INFO [2026-02-26 10:48:28] To cite ggplot2 in publications, please use
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
#> INFO [2026-02-26 10:48:28]   Springer-Verlag New York, 2016.
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28] A BibTeX entry for LaTeX users is
#> INFO [2026-02-26 10:48:28] 
#> INFO [2026-02-26 10:48:28]   @Book{,
#> INFO [2026-02-26 10:48:28]     author = {Hadley Wickham},
#> INFO [2026-02-26 10:48:28]     title = {ggplot2: Elegant Graphics for Data Analysis},
#> INFO [2026-02-26 10:48:28]     publisher = {Springer-Verlag New York},
#> INFO [2026-02-26 10:48:28]     year = {2016},
#> INFO [2026-02-26 10:48:28]     isbn = {978-3-319-24277-4},
#> INFO [2026-02-26 10:48:28]     url = {https://ggplot2.tidyverse.org},
#> INFO [2026-02-26 10:48:28]   }
#> INFO [2026-02-26 10:48:28] Quality metrics were computed as per guidelines in:
#> INFO [2026-02-26 10:48:28] [1] "Broadhurst, David et al. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. Metabolomics : Official journal of the Metabolomic Society vol. 14,6 (2018): 72. doi:10.1007/s11306-018-1367-3"
```
