---
layout: page
title: Cancer analysis
permalink: /4-Cancer.md/
---

Analyzing cancer samples if more challenging. PopSV is able to detect partial signal which would benefit for the discovery of variant in presence of stromal contamination or cell heterogeneity. However we have to be cautious, especially if the genome is highly rearranged.

## More robust normalization

Targeted normalization uses other regions of the genome to normalize one bin. If by chance a majority of these supporting bins are located within copy number aberrations, the normalization might use an incorrect baseline. To avoid this potential over-normalization we implemented a more robust normalization that also take into account the genomic location of the supporting bins and mean coverage in the reference samples.

In practice, this normalization takes a bit more time to run but is recommended when analyzing tumor samples. It can be switched on in `tn.test.sample` function using parameter `aberrant.cases=TRUE`.

## Flagging chromosome aneuploidy

When the a chromosome is not diploid in the tumor, we might want to remove it from the analysis. Indeed, the presence of a large portion of the genome with abnormal coverage will impair the normalization and significance estimation.

In practice we would recommend to start by running PopSV normally. If some samples have a large fraction of affected genome, you could look for aneuploid chromosomes. Once found, these chromosomes could be removed and CNV called again. More subtle variants might be detected.

To flag aneuploid chromosomes in a sample we implemented `aneuploidy.flag` function. It can be run for a specific sample:

```r
aneu.o = aneuploidy.flag("SAMP1827", files.df)
```

Then, `aneu.o$aneu.chrs` contains the names of aneuploid chromosomes and can be added to `call.abnormal.cov` function (see below). In addition, using `plot=TRUE` will produce a graph like this:

*Example of the graph soon...*


## Special parameter when calling abnormal coverage

The fraction of abnormal genome will likely high in a tumor sample compared to normal tissues. For this reason the assumptions used when estimating significance need to be slightly revisited. For example we want to be able to estimate the null distribution even if 30% of the genome is affected.

In practice a few parameters in `call.abnormal.cov` function need to be changed:

+ `min.normal.prop=` represent the minimum proportion of the genome that is assumed to be normal. The default value is `0.9` but we use `0.6` for tumor samples.
+ `aneu.chrs=` can be used to remove aneuploid chromosomes from the test to increase the power in other chromosomes. The vector of aneuploid chromosomes can be computed by `aneuploidy.flag` function (see previous section).
