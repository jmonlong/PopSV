---
layout: page
title: Publications
permalink: /6-Publication.md/
---

## Genome-wide characterization of copy number variants in epilepsy patients

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1002893.svg)](https://doi.org/10.5281/zenodo.1002893)

+ [Manuscript on bioRxiv](https://www.biorxiv.org/content/early/2017/10/06/199224).
+ [GitHub repo](https://github.com/jmonlong/epipopsv) with code to reproduce the figures and numbers.
+ [Data on FigShare](https://figshare.com/s/20dfdedcc4718e465185) includes calls and data necessary to reproduce the analysis.

The CNV calls in the epilepsy patients and controls can download be downloaded from the [FigShare repo](https://figshare.com/s/20dfdedcc4718e465185). The file `cnvs-PopSV-Epilepsy-198affected-301controls-5kb.tsv.gz` (21Mb) is tab-delimited with the following columns:

+ *chr*, *start*, *end*: the position in the reference genome.
+ *project*: either `affected` to denote an epilepsy patient, or `control` for control individuals.
+ *sample*: the sample name.
+ *z*: the Z-score used for calling the CNV.
+ *fc*: the fold change compared to the coverage in the reference samples. <1 for deletions and >1 for duplications.
+ *mean.cov*: the average coverage in the reference samples. 
+ *pv*/*qv*: the P-value and Q-value (or FDR threshold) for the call.
+ *nb.bins.cons*: the number of consecutive 5 Kbp bins in the call.


## Human copy number variants are enriched in regions of low-mappability

+ [Manuscript on bioRxiv](http://www.biorxiv.org/content/early/2015/12/11/034165)
+ [Data and scripts](https://figshare.com/s/8fd3007ebb0fbad09b6d)
+ New version coming up very soon.

The calls from the three datasets (640 genomes total) were merged and formatted into single files to facilitate external analysis/usage of our CNV catalog. The files are available from the [FigShare repository](https://figshare.com/s/8fd3007ebb0fbad09b6d).

- `CNV-PopSV-Twin_CageKid_GoNL-germline.tsv` file contains all the germline CNV calls for the three datasets.
- `CNV-PopSV-CageKid-somatic.tsv` contains the somatic CNV calls from the renal cancer dataset.

Both files are tab-delimited text files with columns:

+ *sample*: the sample name.
+ *chr*, *start*, *end*: the position in the reference genome.
+ *fc*: the fold change compared to the coverage in the reference samples. <1 for deletions and >1 for duplications.
+ *qv*: the Q-value, or FDR threshold, for the call.
+ *cn2.dev*: deviation from the 2 copy number state, i.e. the absolute number of copies of difference between the variant and a diploid state.
+ *cn*: the estimated copy-number state.

Columns *qv* and *cn2.dev* can be used to further select higher confidence calls. For example, one could select calls with small *qv* and/or high *cn2.dev*. Copy number estimation (column *cn*) is only reliable when the the variant spans several bins, hence it's normal to have partial estimates for small calls. Furthermore, copy number estimation in low-coverage regions is more challenging and calls in these regions often exhibits large *cn* (or *fc*).


