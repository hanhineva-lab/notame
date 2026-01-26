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
#> INFO [2026-01-26 11:40:22] Preprocessing and analyses were performed using notame package:
#> INFO [2026-01-26 11:40:22] To cite package ‘notame’ in publications use:
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   Klavus, A.; Kokla, M.; Noerman, S.; Koistinen, V.M.; Tuomainen, M.;
#> INFO [2026-01-26 11:40:22]   Zarei, I.; Meuronen, T.; Hakkinen, M.R.; Rummukainen, S.; Farizah
#> INFO [2026-01-26 11:40:22]   Babu, A.; Sallinen, T.; Karkkainen, O.; Paananen, J.; Broadhurst, D.;
#> INFO [2026-01-26 11:40:22]   Brunius, C.; Hanhineva, K. "notame": Workflow for Non-Targeted LC-MS
#> INFO [2026-01-26 11:40:22]   Metabolic Profiling. Metabolites 2020, 10, 135.
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22] A BibTeX entry for LaTeX users is
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   @Article{,
#> INFO [2026-01-26 11:40:22]     title = {"notame": Workflow for Non-Targeted LC-MS Metabolic Profiling},
#> INFO [2026-01-26 11:40:22]     author = {Anton Klavus and Marietta Kokla and Stefania Noerman and Ville M. Koistinen and Marjo Tuomainen and Iman Zarei and Topi Meuronen and Merja R. Hakkinen and Soile Rummukainen and Ambrin {Farizah Babu} and Taisa Sallinen and Olli Karkkainen and Jussi Paananen and David Broadhurst and Carl Brunius and Kati Hanhineva},
#> INFO [2026-01-26 11:40:22]     year = {2020},
#> INFO [2026-01-26 11:40:22]     journal = {Metabolites},
#> INFO [2026-01-26 11:40:22]     volume = {10},
#> INFO [2026-01-26 11:40:22]     number = {135},
#> INFO [2026-01-26 11:40:22]     url = {https://doi.org/10.3390/metabo10040135},
#> INFO [2026-01-26 11:40:22]     doi = {10.3390/metabo10040135},
#> INFO [2026-01-26 11:40:22]   }
#> INFO [2026-01-26 11:40:22] The primary data structure in notame is SummarizedExperiment:
#> INFO [2026-01-26 11:40:22] To cite package ‘SummarizedExperiment’ in publications use:
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   Morgan M, Obenchain V, Hester J, Pagès H (2025).
#> INFO [2026-01-26 11:40:22]   _SummarizedExperiment: A container (S4 class) for matrix-like
#> INFO [2026-01-26 11:40:22]   assays_. doi:10.18129/B9.bioc.SummarizedExperiment
#> INFO [2026-01-26 11:40:22]   <https://doi.org/10.18129/B9.bioc.SummarizedExperiment>. R package
#> INFO [2026-01-26 11:40:22]   version 1.41.0,
#> INFO [2026-01-26 11:40:22]   <https://bioconductor.org/packages/SummarizedExperiment>.
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22] A BibTeX entry for LaTeX users is
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   @Manual{,
#> INFO [2026-01-26 11:40:22]     title = {SummarizedExperiment: A container (S4 class) for matrix-like assays},
#> INFO [2026-01-26 11:40:22]     author = {Martin Morgan and Valerie Obenchain and Jim Hester and Hervé Pagès},
#> INFO [2026-01-26 11:40:22]     year = {2025},
#> INFO [2026-01-26 11:40:22]     note = {R package version 1.41.0},
#> INFO [2026-01-26 11:40:22]     url = {https://bioconductor.org/packages/SummarizedExperiment},
#> INFO [2026-01-26 11:40:22]     doi = {10.18129/B9.bioc.SummarizedExperiment},
#> INFO [2026-01-26 11:40:22]   }
#> INFO [2026-01-26 11:40:22] visualizations in notame are built with ggplot2:
#> INFO [2026-01-26 11:40:22] To cite ggplot2 in publications, please use
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
#> INFO [2026-01-26 11:40:22]   Springer-Verlag New York, 2016.
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22] A BibTeX entry for LaTeX users is
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   @Book{,
#> INFO [2026-01-26 11:40:22]     author = {Hadley Wickham},
#> INFO [2026-01-26 11:40:22]     title = {ggplot2: Elegant Graphics for Data Analysis},
#> INFO [2026-01-26 11:40:22]     publisher = {Springer-Verlag New York},
#> INFO [2026-01-26 11:40:22]     year = {2016},
#> INFO [2026-01-26 11:40:22]     isbn = {978-3-319-24277-4},
#> INFO [2026-01-26 11:40:22]     url = {https://ggplot2.tidyverse.org},
#> INFO [2026-01-26 11:40:22]   }
data(toy_notame_set)
ex_set <- flag_quality(toy_notame_set)
#> INFO [2026-01-26 11:40:22] 
#> 92% of features flagged for low quality
# Broadhurst et al.(2018) added to citations
citations()
#> INFO [2026-01-26 11:40:22] Preprocessing and analyses were performed using notame package:
#> INFO [2026-01-26 11:40:22] To cite package ‘notame’ in publications use:
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   Klavus, A.; Kokla, M.; Noerman, S.; Koistinen, V.M.; Tuomainen, M.;
#> INFO [2026-01-26 11:40:22]   Zarei, I.; Meuronen, T.; Hakkinen, M.R.; Rummukainen, S.; Farizah
#> INFO [2026-01-26 11:40:22]   Babu, A.; Sallinen, T.; Karkkainen, O.; Paananen, J.; Broadhurst, D.;
#> INFO [2026-01-26 11:40:22]   Brunius, C.; Hanhineva, K. "notame": Workflow for Non-Targeted LC-MS
#> INFO [2026-01-26 11:40:22]   Metabolic Profiling. Metabolites 2020, 10, 135.
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22] A BibTeX entry for LaTeX users is
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   @Article{,
#> INFO [2026-01-26 11:40:22]     title = {"notame": Workflow for Non-Targeted LC-MS Metabolic Profiling},
#> INFO [2026-01-26 11:40:22]     author = {Anton Klavus and Marietta Kokla and Stefania Noerman and Ville M. Koistinen and Marjo Tuomainen and Iman Zarei and Topi Meuronen and Merja R. Hakkinen and Soile Rummukainen and Ambrin {Farizah Babu} and Taisa Sallinen and Olli Karkkainen and Jussi Paananen and David Broadhurst and Carl Brunius and Kati Hanhineva},
#> INFO [2026-01-26 11:40:22]     year = {2020},
#> INFO [2026-01-26 11:40:22]     journal = {Metabolites},
#> INFO [2026-01-26 11:40:22]     volume = {10},
#> INFO [2026-01-26 11:40:22]     number = {135},
#> INFO [2026-01-26 11:40:22]     url = {https://doi.org/10.3390/metabo10040135},
#> INFO [2026-01-26 11:40:22]     doi = {10.3390/metabo10040135},
#> INFO [2026-01-26 11:40:22]   }
#> INFO [2026-01-26 11:40:22] The primary data structure in notame is SummarizedExperiment:
#> INFO [2026-01-26 11:40:22] To cite package ‘SummarizedExperiment’ in publications use:
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   Morgan M, Obenchain V, Hester J, Pagès H (2025).
#> INFO [2026-01-26 11:40:22]   _SummarizedExperiment: A container (S4 class) for matrix-like
#> INFO [2026-01-26 11:40:22]   assays_. doi:10.18129/B9.bioc.SummarizedExperiment
#> INFO [2026-01-26 11:40:22]   <https://doi.org/10.18129/B9.bioc.SummarizedExperiment>. R package
#> INFO [2026-01-26 11:40:22]   version 1.41.0,
#> INFO [2026-01-26 11:40:22]   <https://bioconductor.org/packages/SummarizedExperiment>.
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22] A BibTeX entry for LaTeX users is
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   @Manual{,
#> INFO [2026-01-26 11:40:22]     title = {SummarizedExperiment: A container (S4 class) for matrix-like assays},
#> INFO [2026-01-26 11:40:22]     author = {Martin Morgan and Valerie Obenchain and Jim Hester and Hervé Pagès},
#> INFO [2026-01-26 11:40:22]     year = {2025},
#> INFO [2026-01-26 11:40:22]     note = {R package version 1.41.0},
#> INFO [2026-01-26 11:40:22]     url = {https://bioconductor.org/packages/SummarizedExperiment},
#> INFO [2026-01-26 11:40:22]     doi = {10.18129/B9.bioc.SummarizedExperiment},
#> INFO [2026-01-26 11:40:22]   }
#> INFO [2026-01-26 11:40:22] visualizations in notame are built with ggplot2:
#> INFO [2026-01-26 11:40:22] To cite ggplot2 in publications, please use
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
#> INFO [2026-01-26 11:40:22]   Springer-Verlag New York, 2016.
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22] A BibTeX entry for LaTeX users is
#> INFO [2026-01-26 11:40:22] 
#> INFO [2026-01-26 11:40:22]   @Book{,
#> INFO [2026-01-26 11:40:22]     author = {Hadley Wickham},
#> INFO [2026-01-26 11:40:22]     title = {ggplot2: Elegant Graphics for Data Analysis},
#> INFO [2026-01-26 11:40:22]     publisher = {Springer-Verlag New York},
#> INFO [2026-01-26 11:40:22]     year = {2016},
#> INFO [2026-01-26 11:40:22]     isbn = {978-3-319-24277-4},
#> INFO [2026-01-26 11:40:22]     url = {https://ggplot2.tidyverse.org},
#> INFO [2026-01-26 11:40:22]   }
#> INFO [2026-01-26 11:40:22] Quality metrics were computed as per guidelines in:
#> INFO [2026-01-26 11:40:22] [1] "Broadhurst, David et al. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. Metabolomics : Official journal of the Metabolomic Society vol. 14,6 (2018): 72. doi:10.1007/s11306-018-1367-3"
```
