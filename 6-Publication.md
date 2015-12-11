---
layout: page
title: Publication
permalink: /6-Publication.md/
---

+ [Manuscript on bioRxiv](http://www.biorxiv.org/content/early/2015/12/11/034165)
+ [Data and scripts](https://figshare.com/s/ba79730bb87a1322480d)
+ [PopSV v1.0](https://github.com/jmonlong/PopSV/releases/tag/v1.0) was used in this paper.

# How to reproduce the results from the paper ?

The aligned reads are been deposited by the different projects at different locations:

+ The twins study at the European Nucleotide Archive under [ENA PRJEB8308](https://www.ebi.ac.uk/ena/data/view/PRJEB8308).
+ CageKid (renal cell carcinoma) at the European Genome-phenome Archive under [EGAS00001000083](https://www.ebi.ac.uk/ega/studies/EGAS00001000083).
+ GoNL data can be requested on [their website](http://www.nlgenome.nl/).

All the scripts (mostly R) used to run PopSV on these data are available [there](https://figshare.com/s/ba79730bb87a1322480d). Moreover we uploaded temporary files and scripts to make a re-analysis the easiest possible. Hopefully, anyone can download it, unpack and rerun the different modules to re-compute the graphs, tables and numbers described in the paper.

## Structure of the data/script

*PopSV-masterAnalysis.R* uses functions from *scriptsForPaper.R* to run all the analysis on the different datasets. It creates the graphs and tables used in the manuscript.

The necessary data, as well as the scripts to run the CNV callers, have been packed in *dataRuns.tar*. Run `sh unpackPopSVruns.sh` to unpack it.

*More soon...*


