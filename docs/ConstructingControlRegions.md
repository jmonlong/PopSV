Constructing control regions for enrichment analysis
====================================================

Introduction
------------

When investigating enrichment in specific genomic features, it is common to compare the regions of interest with some controls regions. **A simple approach** to construct control regions could be to **randomly select regions across the genome**. Additionally, it is important for the control regions to have the **same size distribution**.

However, a random distribution across the genome is usually **not realistic**. Likely, regions in the genome were not tested because inaccessible or not included in the analysis. Moreover, **you might want to control for some patterns and look for more**.

For example, we first observed enrichment of CNVs in low-mappability regions. We then wanted to test additional enrichment in different repeat classes. Because repeats are enriched in low-mappability regions, repeats will likely be seen enriched in CNVs. We want to avoid spurious correlation and control for the low-mappability enrichment. By constructing control regions with the same low-mappability enrichment we can can test additional enrichment in the different repeat classes without being biased by the relation between low-mappability regions and repeats.

Constructing control regions with *PopSV* package
-------------------------------------------------

First we load the package and retrieve some annotations to play with.

``` {.r}
library(PopSV)
library(AnnotationHub)
ah = AnnotationHub()
```

    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |                                                                 |   1%
      |                                                                       
      |=                                                                |   1%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |==                                                               |   4%
      |                                                                       
      |===                                                              |   4%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   5%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |====                                                             |   7%
      |                                                                       
      |=====                                                            |   7%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |   8%
      |                                                                       
      |======                                                           |   9%
      |                                                                       
      |======                                                           |  10%
      |                                                                       
      |=======                                                          |  10%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |=======                                                          |  12%
      |                                                                       
      |========                                                         |  12%
      |                                                                       
      |========                                                         |  13%
      |                                                                       
      |=========                                                        |  13%
      |                                                                       
      |=========                                                        |  14%
      |                                                                       
      |=========                                                        |  15%
      |                                                                       
      |==========                                                       |  15%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |===========                                                      |  16%
      |                                                                       
      |===========                                                      |  17%
      |                                                                       
      |===========                                                      |  18%
      |                                                                       
      |============                                                     |  18%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |=============                                                    |  19%
      |                                                                       
      |=============                                                    |  20%
      |                                                                       
      |=============                                                    |  21%
      |                                                                       
      |==============                                                   |  21%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |===============                                                  |  22%
      |                                                                       
      |===============                                                  |  23%
      |                                                                       
      |===============                                                  |  24%
      |                                                                       
      |================                                                 |  24%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |=================                                                |  25%
      |                                                                       
      |=================                                                |  26%
      |                                                                       
      |=================                                                |  27%
      |                                                                       
      |==================                                               |  27%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |===================                                              |  28%
      |                                                                       
      |===================                                              |  29%
      |                                                                       
      |===================                                              |  30%
      |                                                                       
      |====================                                             |  30%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |====================                                             |  32%
      |                                                                       
      |=====================                                            |  32%
      |                                                                       
      |=====================                                            |  33%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |======================                                           |  35%
      |                                                                       
      |=======================                                          |  35%
      |                                                                       
      |=======================                                          |  36%
      |                                                                       
      |========================                                         |  36%
      |                                                                       
      |========================                                         |  37%
      |                                                                       
      |========================                                         |  38%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |=========================                                        |  39%
      |                                                                       
      |==========================                                       |  39%
      |                                                                       
      |==========================                                       |  40%
      |                                                                       
      |==========================                                       |  41%
      |                                                                       
      |===========================                                      |  41%
      |                                                                       
      |===========================                                      |  42%
      |                                                                       
      |============================                                     |  42%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |============================                                     |  44%
      |                                                                       
      |=============================                                    |  44%
      |                                                                       
      |=============================                                    |  45%
      |                                                                       
      |==============================                                   |  45%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |==============================                                   |  47%
      |                                                                       
      |===============================                                  |  47%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |================================                                 |  50%
      |                                                                       
      |=================================                                |  50%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |=================================                                |  52%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |==================================                               |  53%
      |                                                                       
      |===================================                              |  53%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |===================================                              |  55%
      |                                                                       
      |====================================                             |  55%
      |                                                                       
      |====================================                             |  56%
      |                                                                       
      |=====================================                            |  56%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |=====================================                            |  58%
      |                                                                       
      |======================================                           |  58%
      |                                                                       
      |======================================                           |  59%
      |                                                                       
      |=======================================                          |  59%
      |                                                                       
      |=======================================                          |  60%
      |                                                                       
      |=======================================                          |  61%
      |                                                                       
      |========================================                         |  61%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |=========================================                        |  62%
      |                                                                       
      |=========================================                        |  63%
      |                                                                       
      |=========================================                        |  64%
      |                                                                       
      |==========================================                       |  64%
      |                                                                       
      |==========================================                       |  65%
      |                                                                       
      |===========================================                      |  65%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |============================================                     |  67%
      |                                                                       
      |============================================                     |  68%
      |                                                                       
      |=============================================                    |  68%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |=============================================                    |  70%
      |                                                                       
      |==============================================                   |  70%
      |                                                                       
      |==============================================                   |  71%
      |                                                                       
      |==============================================                   |  72%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |===============================================                  |  73%
      |                                                                       
      |================================================                 |  73%
      |                                                                       
      |================================================                 |  74%
      |                                                                       
      |================================================                 |  75%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |=================================================                |  76%
      |                                                                       
      |==================================================               |  76%
      |                                                                       
      |==================================================               |  77%
      |                                                                       
      |==================================================               |  78%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |===================================================              |  79%
      |                                                                       
      |====================================================             |  79%
      |                                                                       
      |====================================================             |  80%
      |                                                                       
      |====================================================             |  81%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |=====================================================            |  82%
      |                                                                       
      |======================================================           |  82%
      |                                                                       
      |======================================================           |  83%
      |                                                                       
      |======================================================           |  84%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |=======================================================          |  85%
      |                                                                       
      |========================================================         |  85%
      |                                                                       
      |========================================================         |  86%
      |                                                                       
      |========================================================         |  87%
      |                                                                       
      |=========================================================        |  87%
      |                                                                       
      |=========================================================        |  88%
      |                                                                       
      |==========================================================       |  88%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |==========================================================       |  90%
      |                                                                       
      |===========================================================      |  90%
      |                                                                       
      |===========================================================      |  91%
      |                                                                       
      |===========================================================      |  92%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |============================================================     |  93%
      |                                                                       
      |=============================================================    |  93%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |=============================================================    |  95%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |==============================================================   |  96%
      |                                                                       
      |===============================================================  |  96%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |===============================================================  |  98%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |================================================================ |  99%
      |                                                                       
      |=================================================================|  99%
      |                                                                       
      |=================================================================| 100%

``` {.r}
genes = ah[["AH49010"]] ## Genes
dgv = ah[["AH5120"]] ## SVs from DGV
dgv = dgv[sample.int(length(dgv), 1e4)] ## Reduce to 10K random SVs
```

We imported a gene annotation and 10 thousands SVs from DGV. If we want to construct control regions that fit the SV size and overlap with genes, we run:

``` {.r}
dgv.cont = draw.controls(dgv, list(gene=genes), chr.prefix="chr")
```

Now let's verify that the size distribution is the same.

``` {.r}
library(ggplot2)
size.df = rbind(data.frame(reg="dgv", size=width(dgv)),
                 data.frame(reg="control", size=width(dgv.cont)))
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge")
```

![](ConstructingControlRegions_files/figure-markdown_github/unnamed-chunk-3-1.png)<!-- -->

``` {.r}
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge") + scale_x_log10()
```

![](ConstructingControlRegions_files/figure-markdown_github/unnamed-chunk-3-2.png)<!-- -->

And that the input and output regions overlap genes similarly.

``` {.r}
mean(overlapsAny(dgv, genes))
```

    ## [1] 0.3964

``` {.r}
mean(overlapsAny(dgv.cont, genes))
```

    ## [1] 0.3979

`draw.controls` functions can **accept any number of genomic features to control**. Let's import two additional genomic annotation that we would like to control for our enrichment analysis: assembly gaps and segmental duplications.

``` {.r}
gap = ah[["AH6444"]]
segdups = ah[["AH5121"]]
dgv.cont2 = draw.controls(dgv, list(gene=genes, gap=gap, sd=segdups), chr.prefix="chr")
```

Again, the size distribution must be the same:

``` {.r}
## Same size distribution ?
size.df = rbind(data.frame(reg="dgv", size=width(dgv)),
                 data.frame(reg="control", size=width(dgv.cont2)))
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge")
```

![](ConstructingControlRegions_files/figure-markdown_github/unnamed-chunk-6-1.png)<!-- -->

``` {.r}
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge") + scale_x_log10()
```

![](ConstructingControlRegions_files/figure-markdown_github/unnamed-chunk-6-2.png)<!-- -->

And the overlap with the three different genomic annotations similar.

``` {.r}
## Same overlap with features ?
mean(overlapsAny(dgv, genes))
```

    ## [1] 0.3964

``` {.r}
mean(overlapsAny(dgv.cont2, genes))
```

    ## [1] 0.4008

``` {.r}
mean(overlapsAny(dgv, gap))
```

    ## [1] 0.0063

``` {.r}
mean(overlapsAny(dgv.cont2, gap))
```

    ## [1] 0.0068

``` {.r}
mean(overlapsAny(dgv, segdups))
```

    ## [1] 0.2114

``` {.r}
mean(overlapsAny(dgv.cont2, segdups))
```

    ## [1] 0.21

If we had used the first set of control regions (only genes overlap control) the gap and segmental duplication overlap proportions wouldn't match.

``` {.r}
mean(overlapsAny(dgv.cont, segdups))
```

    ## [1] 0.0924

``` {.r}
mean(overlapsAny(dgv.cont, gap))
```

    ## [1] 0.1019

Controlling for the distance to a feature
-----------------------------------------

In addition to controlling for the overlap to a set of feature, `draw.controls` can also **control for the distance to one feature**. Although we can control for overlap to several feature the distance control is more complex to multiplex and for now we can only control for distance to one feature.

``` {.r}
dgv.cont3 = draw.controls(dgv, list(gene=genes, sd=segdups), chr.prefix="chr", dist.gr=gap)
```

Again, the size distribution must be the same:

``` {.r}
## Same size distribution ?
size.df = rbind(data.frame(reg="dgv", size=width(dgv)),
                 data.frame(reg="control", size=width(dgv.cont3)))
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge")
```

![](ConstructingControlRegions_files/figure-markdown_github/unnamed-chunk-10-1.png)<!-- -->

``` {.r}
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge") + scale_x_log10()
```

![](ConstructingControlRegions_files/figure-markdown_github/unnamed-chunk-10-2.png)<!-- -->

And the overlap with the three different genomic annotations similar.

``` {.r}
## Same overlap with features ?
mean(overlapsAny(dgv, genes))
```

    ## [1] 0.3964

``` {.r}
mean(overlapsAny(dgv.cont3, genes))
```

    ## [1] 0.4016

``` {.r}
mean(overlapsAny(dgv, segdups))
```

    ## [1] 0.2114

``` {.r}
mean(overlapsAny(dgv.cont3, segdups))
```

    ## [1] 0.2091

Finally let's check that we now control for the *distance* to a gap.

``` {.r}
## Same distance to gap ?
dist.df = rbind(data.frame(reg="dgv", dist=as.data.frame(distanceToNearest(dgv, gap))$distance),
    data.frame(reg="simple control", dist=as.data.frame(distanceToNearest(dgv.cont2, gap))$distance),
    data.frame(reg="distance control", dist=as.data.frame(distanceToNearest(dgv.cont3, gap))$distance))
ggplot(dist.df, aes(x=dist, fill=reg)) + geom_histogram(position="dodge") + xlab("distance to a gap")
```

![](ConstructingControlRegions_files/figure-markdown_github/unnamed-chunk-12-1.png)<!-- -->

``` {.r}
ggplot(dist.df, aes(x=dist, colour=reg)) + stat_ecdf() + ylab("cumulative proportion of regions") + xlab("distance to a gap")
```

![](ConstructingControlRegions_files/figure-markdown_github/unnamed-chunk-12-2.png)<!-- -->

R session
---------

``` {.r}
sessionInfo()
```

    ## R version 3.2.3 (2015-12-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 14.04.4 LTS
    ## 
    ## locale:
    ##  [1] LC_CTYPE=fr_CA.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=fr_CA.UTF-8        LC_COLLATE=fr_CA.UTF-8    
    ##  [5] LC_MONETARY=fr_CA.UTF-8    LC_MESSAGES=fr_CA.UTF-8   
    ##  [7] LC_PAPER=fr_CA.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=fr_CA.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_2.1.0        GenomicRanges_1.22.4 GenomeInfoDb_1.6.3  
    ## [4] IRanges_2.4.8        S4Vectors_0.8.11     AnnotationHub_2.2.3 
    ## [7] BiocGenerics_0.16.1  PopSV_1.0            rmarkdown_0.9.5     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.3                       plyr_1.8.3                       
    ##  [3] futile.logger_1.4.1               BiocInstaller_1.20.1             
    ##  [5] formatR_1.2.1                     XVector_0.10.0                   
    ##  [7] futile.options_1.0.0              bitops_1.0-6                     
    ##  [9] tools_3.2.3                       zlibbioc_1.16.0                  
    ## [11] digest_0.6.9                      gtable_0.2.0                     
    ## [13] RSQLite_1.0.0                     evaluate_0.8                     
    ## [15] BSgenome_1.38.0                   shiny_0.13.1                     
    ## [17] DBI_0.3.1                         curl_0.9.6                       
    ## [19] yaml_2.1.13                       rtracklayer_1.30.2               
    ## [21] httr_1.1.0                        stringr_1.0.0                    
    ## [23] knitr_1.12.3                      Biostrings_2.38.4                
    ## [25] grid_3.2.3                        data.table_1.9.6                 
    ## [27] Biobase_2.30.0                    R6_2.1.2                         
    ## [29] BSgenome.Hsapiens.UCSC.hg19_1.4.0 AnnotationDbi_1.32.3             
    ## [31] XML_3.98-1.4                      BiocParallel_1.4.3               
    ## [33] lambda.r_1.1.7                    magrittr_1.5                     
    ## [35] scales_0.4.0                      GenomicAlignments_1.6.3          
    ## [37] Rsamtools_1.22.0                  htmltools_0.3                    
    ## [39] SummarizedExperiment_1.0.2        colorspace_1.2-6                 
    ## [41] mime_0.4                          interactiveDisplayBase_1.8.0     
    ## [43] xtable_1.8-2                      httpuv_1.3.3                     
    ## [45] labeling_0.3                      stringi_1.0-1                    
    ## [47] munsell_0.4.3                     RCurl_1.95-4.8                   
    ## [49] chron_2.3-47
