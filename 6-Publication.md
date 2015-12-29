---
layout: page
title: Publication
permalink: /6-Publication.md/
---

+ [Manuscript on bioRxiv](http://www.biorxiv.org/content/early/2015/12/11/034165)
+ [Data and scripts](https://figshare.com/s/ba79730bb87a1322480d)
+ [PopSV v1.0](https://github.com/jmonlong/PopSV/releases/tag/v1.0) was used in this paper.

All the scripts (mostly R) used to run PopSV on these data are available [there](https://figshare.com/s/ba79730bb87a1322480d). Moreover we uploaded temporary files and scripts to make a re-analysis the easiest possible. Hopefully, anyone can download it, unpack and rerun the different modules to re-compute the graphs, tables and numbers described in the paper.

Here are some details on the [data in FigShare](https://figshare.com/s/ba79730bb87a1322480d).

## Summary of the content

To download PopSV calls in the three datasets, go directly to the last section *CNV catalog*.

To review or reproduce the graphs and numbers from the paper, extract and set up the data (next section), then jump directly to section *Running the analysis*.

To review or reproduce all the runs and analysis from scratch, read the sections in order.

## Getting the data and setting up the structure

The aligned reads (BAM files) are not hosted here but were deposited in different repositories:

+ The twins study at the European Nucleotide Archive under [ENA PRJEB8308](https://www.ebi.ac.uk/ena/data/view/PRJEB8308).
+ CageKid (renal cell carcinoma) at the European Genome-phenome Archive under [EGAS00001000083](https://www.ebi.ac.uk/ega/studies/EGAS00001000083).
+ GoNL data can be requested on [their website](http://www.nlgenome.nl/).

However the bin counts, i.e. the number of reads mapping in each bin of a binned genome, is provided here for the Twin study and our renal cancer dataset. That way, we can reproduce the entire analysis, skipping only the first step: counting the reads from the BAM files.

The bin counts for the different datasets and bin sizes have been archived into `dataRuns.tar` for convenience. We provide an unpacking script to facilitate setting up the data. After downloading the data, run `sh unpackPopSVruns.sh` to unpack the data.

Once unpacked the data and scripts are organized as the following:

+ `Twins` and `CageKid` folders contain raw data and calls on each of these datasets.
+ `annotation` folder contains a script to download and format relevant genomic annotations.
+ `installArchive` folder contains the install files for *FREEC* and *PopSV*.
+ `PopSV-masterAnalysis.R` and `scriptsForPaper.R` are the scripts that produce all the graphs and numbers used in the paper.

### Setting up the annotation files

In the `annotation` folder, run:

+ `sh downloads.sh` to download annotations (mainly from UCSC).
+ `Rscript dataFormat.R` to format them for easy use later in the analysis. Of note, this step requires packages *GenomicRanges* and
*dplyr*.

### Installing the softwares

#### PopSV
The version of PopSV used in our paper is available in folder `installArchive` and can be installed as a normal R package using:
```
R CMD INSTALL PopSV_1.0.tar.gz
```

#### FREEC
To install FREEC, follow instruction on their [website](http://bioinfo-out.curie.fr/projects/freec/). The version we used can be found in `installArchive` folder.

#### cn.MOPS
We installed cn.MOPS through [Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/cn.mops.html). We used version 1.12.0.


## Running PopSV on the two datasets

**This part is optional.** In order to reproduce the graphs and numbers form the paper, there is NO need to rerun PopSV again. If you are only interested in reproducing the graphs starting from PopSV output, the necessary files are also deposited here and you can go directly to the *Running the analysis* section further down. The scripts described now are present for transparency purpose or if someone wants to rerun the full analysis from scratch.

Each PopSV run, i.e. for a dataset and a specific bin size, has been run following the script called `run-PopSV-<DATASET>-<BINSIZE>.R`. For example, the script used to run the 5 Kbp bin analysis on the Twins study is called `run-PopSV-twin-5kbp.R` and is located in `Twins/5kbp/` folder.

These scripts runs PopSV step by step and use *BatchJobs* package to send the computations to a computing cluster. More details on the different steps and usage can be found [there](http://jmonlong.github.io/PopSV//1-BasicWorkflow.md/) and [there](http://jmonlong.github.io/PopSV//2-ClusterManagement.md/).


## Running FREEC and cn.MOPS

**This part is optional.** Same here, you don't have to rerun these methods if you just want to analyze their results. Jump directly to the next section if you are not interested in running the full analysis from scratch.

Located in `Twins/otherMethods/FREEC/` folder, the two scripts `run-freec-twins-500bp.sh` and `run-freec-twins-5kbp.sh` show how CNVs were called using FREEC. Similarly, the scripts for CageKid dataset are located in `CageKid/otherMethods/FREEC/`.

Located in `Twins/otherMethods/cn.MOPS/` folder, the  `run-cnmops-twins.sh` shows how CNVs were called using cn.MOPS. Similarly, the scripts for CageKid dataset are located in `CageKid/otherMethods/cn.MOPS/`.


## Running the analysis

We also deposited the files necessary to rerun the analysis and produce graphs and numbers from the paper. These files are located at different locations and could be reproduced using the instructions from the two previous sections. For example:

+ `cnvs-twin-FDR001-mergeStitch-thSdest.RData` and `otherMethods-500bp-5kbp.RData` in `Twins/` contains the final calls from PopSV and the merged calls from FREEC and cn.MOPS.
+ `cnvs-twin-5kbp-FDR001-mergeStitch-thSdest.RData` in `Twins/5kbp/` contains the calls from PopSV for the 5 Kbp bin run.
+ `ref-msd.tsv` in `Twins/5kbp/` contains the average coverage in the reference samples across the genome.
+ Same for CageKid in `CageKid/` folder.

`PopSV-masterAnalysis.R` and `scriptsForPaper.R` can be used to rerun the analysis described in the paper.

First, `scriptsForPaper.R` contains functions that each performs one specific analysis. For example, function *biasWGS* measures the bias present in WGS data, function *concordanceTwins* compare the calls between twins and computes quality measures. These functions are generic enough so that we can run them on different datasets or different bin size runs, without duplicating code. And that's exactly what `PopSV-masterAnalysis.R` script does: it loads `scriptsForPaper.R` and calls the functions in the different datasets and bin sizes.

Before running anything, you should check that all the necessary packages are installed. Try running `source("scriptsForPaper.R")` in R. If there is no error, you have all the necessary packages. Else, install the missing package (see error message) using either `install.packages("theMissingPackage")` or [Bioconductor install function](http://bioconductor.org/install/#install-bioconductor-packages).

In `PopSV-masterAnalysis.R` script, a few functions are first defined to help automating the analysis. Then the analysis is grouped by dataset, and then by bin size. At the very bottom, you can find the analysis involving all datsets, e.g. summary tables. A few notes:

+ `NB.CORES = 3` can be changed if you want to use more/less computing nodes.
+ Individual functions can be run interactively, i.e. no need to run the entire script every time.
+ If you run only a particular function don't forget to first run the few lines on top of its paragraph (up to the `attach(..)` function). These lines define the data to be used by the function. Also don't forget the first part of the script which loads `scriptsForPaper.R`, defines useful functions and data paths.

Last note, you might notice that the function are not directly called but instead called through `bjify(...)`. The short explanation for that is the poor memory management in R. I use *BatchJobs* to open R sub-processes, run the function, get the output and close the sub-process. This hack avoids having a memory explosion when running several functions one after the other. As a side note, it shouldn't be like this, one of the motivation to use separate functions for separate analysis was exactly to avoid memory explosion, but apparently the function environments are not cleaned properly and this (ugly) hack was necessary. R can be disappointing sometimes, but *BatchJobs* saved the day. Here, there is no real need to configure *BatchJobs* as it runs locally. Eventually, create a `.BatchJobs.R` file with `cluster.functions <- makeClusterFunctionsLocal()` in it, and it should work.


## CNV catalog

The calls from the three datasets were merged and formatted into single files to facilitate external analysis/usage of our CNV catalog. `CNV-PopSV-Twin_CageKid_GoNL-germline.tsv` file contains all the germline CNV calls for the three datasets. `CNV-PopSV-CageKid-somatic.tsv` contains the somatic CNV calls from the renal cancer dataset.

Both files are tab-delimited text files with columns:

+ *sample*: the sample name.
+ *chr*, *start*, *end*: the position in the reference genome.
+ *fc*: the fold change compared to the coverage in the reference samples. <1 for deletions and >1 for duplications.
+ *qv*: the Q-value, or FDR threshold, for the call.
+ *cn2.dev*: deviation from the 2 copy number state, i.e. the absolute number of copies of difference between the variant and a diploid state.
+ *cn*: the estimated copy-number state.

Columns *qv* and *cn2.dev* can be used to further select higher confidence calls. For example, one could select calls with small *qv* and/or high *cn2.dev*. Copy number estimation (column *cn*) is only reliable when the the variant spans several bins, hence it's normal to have partial estimates for small calls. Furthermore, copy number estimation in low-coverage regions is more challenging and calls in these regions often exhibits large *cn* (or *fc*).

