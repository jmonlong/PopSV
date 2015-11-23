---
layout: page
title: Frequently Asked Questions
permalink: /faq/
---

*Soon...*

### What bin size should I use ?

Depends on your sequencing depth...

### Can I run PopSV on my laptop ?

Yes but you shouldn't...

### Do I need permission to install PopSV on the cluster ?

No, you can install the package in you home folder if you want. For example to install it in `~/R` , use:

```r
devtools::install_github("jmonlong/PopSV", ref="forPaperXY", args="-l ~/R")
```

Then you should load the package with:

```r
library(PopSV, lib.loc="~/R")
```

### Are `FDR.th=` parameter and `qv` column the same thing ?

Yes and no...


### X and Y chromosomes are not tested by default, why ? I want them !

By default PopSV analyzes all samples together. In practice males and females are jointly analyzed so we only analyze chromosomes where they are expected to have similar copies, i.e. the autosomes.

However, if needed, the samples can be split into males and females, and chromosomes X and Y analyzed. See [this script](https://github.com/jmonlong/PopSV/blob/forPaper/scripts/run-PopSV-batchjobs-XY-automatedPipeline.R) for an example of what the pipeline would look like. There three analysis are run:

1. All samples on chr 1 to 22.
2. Females on chr X.
3. Males on chr X and Y.
