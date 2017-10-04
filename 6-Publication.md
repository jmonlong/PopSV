---
layout: page
title: Publication
permalink: /6-Publication.md/
---

## Human copy number variants are enriched in regions of low-mappability

+ [Manuscript on bioRxiv](http://www.biorxiv.org/content/early/2015/12/11/034165)
+ [Data and scripts](https://figshare.com/s/ba79730bb87a1322480d)
+ New version coming up very soon.

The calls from the three datasets (640 genomes total) were merged and formatted into single files to facilitate external analysis/usage of our CNV catalog. The files are available from the [FigShare repository](https://figshare.com/s/ba79730bb87a1322480d).

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

