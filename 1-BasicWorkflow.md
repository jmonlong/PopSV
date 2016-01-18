---
layout: page
title: Analysis steps
permalink: /1-BasicWorkflow.md/
---

Examples of an analysis, for local computation or computing cluster usage, can be found on the [`scripts` folder of the GitHub repository](https://github.com/jmonlong/PopSV/tree/master/scripts). For more information on a specific function, see the manual or access the documentation through `?the.function.name`.

## Input files
The analysis can start directly from the BAM files. Each BAM file needs to be **sorted** and **indexed** (see [samtools](http://www.htslib.org/)).

A tabular separated values (*tsv*) file with the name of the sample (in column named *sample*) and the path to the corresponding BAM (in column named *bam*) is imported and  given to the `init.filenames`. This function will create the path and file names for the different files created and used throughout the analysis.

Then, the regions of interest or *bins* have to be defined. `fragment.genome.hp19` automatically fragments hg19 genome into non-overlapping consecutive windows of a specified size. However PopSV can perform with any type of windows. It is still recommended to define non-overlapping windows and, for computation reasons, no more than a total of 10 million of bins. If a custom definition is used it should be a *data.frame* with columns named *chr*, *start* and *end*.

Finally, the GC content of each bin can be computed, for hg19, using function `getGC.hg19`. If another genome is to be used, GC content should be defined as a *data.frame* with columns named *chr*, *start*, *end* and *GCcontent*.

## Counting reads
Reads are counted in each bin to measure coverage. `bin.bam` function will count reads in each bin for a given sample.

Eventually, this can be done externally, e.g. using [bedtools coverage](http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html). The final count file should have these four columns: *chr*, *start*, *end* and *bc* (for *b*in *c*ount).

## GC bias correction
GC bias is corrected using a LOESS model. Using this model, a normalization coefficient is computed for each bin based on its GC content. This step is performed by `correct.GC` function.

## Sample Quality Control
The last "pre-processing" step aims at defining the set of reference samples. These samples will define "normal" coverage. A natural set of reference samples are the controls in case/control studies or normal samples in normal/tumor paired samples designs. The more the better but **reference samples should be homogeneous** to get optimal detection power. Eventually, if all available samples are normal, they can all be used as reference and later tested against themselves. 

`qc.samples.cluster` opens an interactive web-browser application to explore how homogeneous the samples are. The samples are represented and clustered using the first two principal components. You can then decide how many clear clusters are present. If any the different batches could be analyzed separately or some clusters completely removed (outliers). **This step is usually not necessary** because samples analyzed should have been sequenced with similar protocol and technologies. It is mostly a safety check. For this QC, a subset of the bins can be used, e.g. using `quick.count` function. More information on the usage [there]({{ site.baseurl }}3-Visualization.md#data-quality-before-analysis).

`qc.samples` function will join all bin count files and produce some graphs on the set of reference samples. If too many reference samples are available (lucky you), `nb.ref.samples=` parameters can be used. 200 reference samples is usually enough.

## Normalization of the reference samples

This step is extremely important to avoid sample-specific bias being picked up as abnormal coverage. Targeted normalization is implemented in the `tn.norm` function. This normalization works better than more naive/global approaches (in my experience). Other methods might be integrated at some point (don't hesitate to [ask](https://github.com/jmonlong/PopSV/issues)). For example Principal Component removal, quantile normalization, or more naive normalization of the coverage median and variance, might be enough (maybe) if the data quality is very good.

The output of the normalization function is a *data.frame* with the normalized bin counts as well as some metrics on the normalization efficiency.

## Z-score computation for the reference samples

From the normalized bin counts the mean and standard deviation across reference samples is computed for each bin. A Z-score and fold-change estimate are derived for each reference sample. Function `z.comp` performs this step.

## Testing other samples

At this point the rest of the samples can be tested using `tn.test.sample` function. Both normalization and Z-score computation are performed by the function.

## Abnormal coverage calls
To find which bins have abnormally low/high coverage, `call.abnormal.cov` will derive P-value from the Z-score from a particular sample and use False Discovery Rate control. A stitching bin merging strategy is performed (parameters `merge.cons.bins=` and `stitch.dist=`). Other parameters are described in the function documentation, such as the strategy for P-value definition or extra tricks (e.g. [for cancer samples]({{ site.baseurl }}4-Cancer.md)). 

## Visualizing the results
Function `sv.summary.interactive` opens an interactive web-browser application to look at the results. The number of calls across samples, distribution of copy-number estimates or frequency is visualized. You can play with different stringency filters in order to retrieve the ones that gives the best result quality.

More on the different visualizations available in PopSV [here]({{ site.baseurl }}3-Visualization.md).
