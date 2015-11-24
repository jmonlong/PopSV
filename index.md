---
layout: default
title: PopSV
---

# What is PopSV ?

PopSV is a Copy-Number Variation (CNV) detection method from high-throughput sequencing. Abnormal Read-Depth signal is detected by using a population of samples as reference. Thanks to this population view the whole genome can be robustly interrogated, including regions of low mappability. Moreover, any divergence from the reference samples are detected, even if the signal is incomplete, e.g. tumoral aberrations or SV involving repeats.

**Warning: PopSV package is still in active development.**

# Getting started

## Installation

To install the latest development version:

```r
devtools::install_github("jmonlong/PopSV")
```

This command requires [*devtools* package](https://github.com/hadley/devtools) which can be easily installed with :

```r
install.packages("devtools")
```

Some [Bioconductor](http://bioconductor.org/) packages are also necessary and not installed automatically. Running the following command should be sufficient :

```r
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19", "Rsamtools")
```

**R 3.1 or higher** is required.

## Usage

PopSV package is used as any R package, by simply loading the library and using the provided funtions.

```r
library(PopSV)
...
```

## I don't want to think, how can I run it quickly ?

After quickly configuring *BatchJobs* for your computing cluster, you can run the [automated pipeline]({{ site.baseurl }}2-ClusterManagement.md#automated-run):

```r
## Load package and wrapper
library(BatchJobs)
library(PopSV)
source("automatedPipeline.R")
## Set-up files and bins
bam.files = read.table("bams.tsv", as.is=TRUE, header=TRUE)
files.df = init.filenames(bam.files, code="example")
save(files.df, file="files.RData")
bin.size = 1e3
bins.df = fragment.genome.hp19(bin.size)
save(bins.df, file="bins.RData")
## Run PopSV
res.GCcounts = autoGCcounts("files.RData", "bins.RData")
res.df = autoNormTest("files.RData", "bins.RData")
```



# Workflow

![PopSV workflow](public/PopSVworkflow.png)

First the genome is fragmented and reads mapping in each bin are counted for each sample and GC corrected (1). Next, coverage of the sample is normalized (2) and each bin is tested by computing a Z-score (3), estimating p-values (4) and identifying abnormal regions (5).

A quick description of the different analysis steps and their corresponding functions can be found in [this page]({{ site.baseurl }}1-BasicWorkflow.md).


# Running PopSV on computing clusters

PopSV workflow heavily uses [*BatchJobs* package](https://cran.r-project.org/web/packages/BatchJobs/index.html) to send computations to a cluster. It needs some configuration but then it saves a lot of time and the pipeline can be run easily. For more information on how to configure it and how the pipelines are using it go to [this page]({{ site.baseurl }}2-ClusterManagement.md).


# FAQ

Find frequently asked questions [there]({{ site.baseurl }}4-FAQ.md).
