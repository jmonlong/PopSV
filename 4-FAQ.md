---
layout: page
title: Frequently Asked Questions
permalink: /4-FAQ.md/
---


### Can I run PopSV on my laptop ?

Yes, it's possible but not recommended. You can have a look at the [pipeline to run PopSV locally](https://github.com/jmonlong/PopSV/blob/forPaper/scripts/run-PopSV-local.R). However, PopSV has been designed to be easily run in a cluster using *BatchJobs* package. Moreover, you likely have access to a computing cluster to store and analyze several WGS samples. Have a look [there]({{ site.baseurl }}2-ClusterManagement.md#installation-and-configuration) to see how to configure it.

Eventually if you really want to run it locally (and slowly) but want to benefit from the optimization with *BatchJobs* (e.g. the [automated pipeline]({{ site.baseurl }}2-ClusterManagement.md#automated-run)), you can use your laptop/computer as cluster. Then you don't need [3 configuration files]({{ site.baseurl }}2-ClusterManagement.md#installation-and-configuration), just create a `.BatchJobs.R` file in the project folder or `~/` with

```r
cluster.functions <- makeClusterFunctionsMulticore(ncpus=6)
```

With this configuration 6 cores of your computer will be used to send jobs.


### Do I need permission to install PopSV on the cluster ?

No, you can install the package in you home folder if you want. For example to install it in a `~/R` folder, use:

```r
devtools::install_github("jmonlong/PopSV", args="-l ~/R")
```

Then you should load the package with:

```r
library(PopSV, lib.loc="~/R")
```

Eventually, if you are using the [automated pipeline]({{ site.baseurl }}2-ClusterManagement.md#automated-run) you can pass this pass to `lib.loc=` parameter.


### What bin size should I use ?

It depends on your sequencing depth. We usually perform a run with a small bin and one with a large bin. In practice we aim at an average coverage of 200/2000 reads per bin. For ~100 bp reads with a sequencing depth of 40X, we then use 500 bp (200*100/40) and 5 Kbp.

To a certain point, the homogeneity of your samples would also affect what is the best bin size. More heterogeneity, i.e. technical noise, would require slightly larger bins.

### How do I know if my reference samples are correct ?

First, the reference samples should be normal genomes. They will be used as baseline to define what is a normal genome so in most cases they should be healthy individuals.

Second, they need to be homogeneous. `qc.samples.summary` function can be used to interactively cluster the samples using the coverage. It should be done right after counting the reads in all the samples. More details [there]({{ site.baseurl }}1-BasicWorkflow.md#sample-quality-control).


### Are `FDR.th=` parameter and `qv` column the same thing ?

Yes and no. They represent the same measure but the `qv` column is a merged across all the bins of a call, while `FDR.th` in `call.abnormal.cov` function is used to select each bin before merging them into calls.

In summary, it more powerful to play with `FDR.th=` parameter than to filter directly the `qv` column. Better run `call.abnormal.cov` several times than running it once filtering the output.

More details on `call.abnormal.cov` function [there]({{ site.baseurl }}1-BasicWorkflow.md#abnormal-coverage-calls).

### X and Y chromosomes are not tested, why ? I want them !

By default PopSV analyzes all samples together. In practice males and females are jointly analyzed so we only analyze chromosomes where they are expected to have similar copies, i.e. the autosomes.

However, if needed, the samples can be split into males and females, and chromosomes X and Y analyzed. See [this script](https://github.com/jmonlong/PopSV/blob/forPaper/scripts/run-PopSV-XY-batchjobs-automatedPipeline.R) for an example of what the pipeline would look like. There three analysis are run:

1. All samples on chr 1 to 22.
2. Females on chr X.
3. Males on chr X and Y.
