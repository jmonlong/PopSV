PopSV
=====
PopSV is a structural variation (SV) detection method from high-throughput sequencing. 
Abnormal Read-Depth signal is detected by using a population of samples ase reference. Thanks to this population
view the whole genome can be robustly interrogated, including regions of low mappability. Moreover, any divergence from
the reference samples are detected, even if the signal is incomplete, e.g. tumoral aberrations or SV involving repeats.

**Warning: PopSV package is still in active development.**

### Installation
To install the latest development version: `devtools::install_github("jmonlong/PopSV")`. This requires `devtools` package (more information [here](https://github.com/hadley/devtools)) which can be easily installed with `install.packages("devtools")`. 

It requires R 3.1 or higher.

### Analysis steps
Examples of a analysis, for local computation or computing cluster usage, can be found on the `scripts` folder. For more information on a specific function, see the manual or access the documentation through `?the.function.name`.

#### Input files
The analysis can start directly from the BAM files. Each BAM file needs to be **sorted** and **indexed** (see [samtools](http://www.htslib.org/)).

A tabular separated values (*tsv*) file with the name of the sample (in column named *sample*) and the path to the corresponding BAM (in column named *bam*) is imported and  given to the `init.filenames`. This function will create the path and file names for the different files created and used throughout the analysis.

Then, the regions of interest or *bins* have to be defined. `fragment.genome.hp19` automatically fragments hg19 genome into non-overlapping consecutive windows of a specified size. However PopSV can perform with any type of windows. It is still recommended to define non-overlapping windows and, for computation reasons, no more than a total of several million of bins. If a custom definition is used it should follow the BED format with columns named *chr*, *start* and *end*.

Finally, the GC content of each bin can be computed, for hg19, using function `getGC.hg19`. If another genome is to be used, GC content should be define following BedGraph format, i.e with columns named *chr*, *start*, *end* and *GCcontent*.

#### Counting reads
Reads are counted in each bin to measure coverage. `bin.bam` function will count reads in each bin for a given sample.

Eventually, this can be done externally, e.g. using [bedtools coverage](). The final count file should have these four columns: *chr*, *start*, *end* and *bc* (for bin count).

#### GC bias correction
GC bias is corrected using a LOESS model. Using this model, a normalization coefficient is computed for each bin based on its GC content. This step is performed by `correct.GC` function.

#### Sample Quality Control
The last "pre-processing" step aims at defining the set of reference samples. These samples will define "normal" coverage. A natural set of reference samples are the controls in case/control studies or normal samples in normal/tumor paired samples designs. The more the better but these samples should also be homogeneous to get optimal detection power.

`qc.samples` function will join all bin count files and compute some comparative measures to help define an homogeneous set of reference samples. Using the output of `qc.samples`, `qc.samples.summary` open an interactive visualization application, in a web browser, to best assess which threshold to choose. 

#### Normalization and Z-score computation
This is the core of the test where the coverage in each bin is normalized and tested against the reference. Depending on the normalization method used different function can be used: `tn.norm` for Targeted Normalization (recommended), `pca.norm` for Principal Component regression normalization, `medvar.norm` for a more naive normalization of the coverage median and variance.

The output of these functions is a list with:

+ *z*: the Z-scores for each bin and each sample.
+ *cn.coeff*: an estimated copy number for each bin and each sample.
+ *bc.norm*: the normalized read count for each bin and each sample.
+ *norm.stats*: some metrics on the efficiency of the normalization.

#### Abnormal coverage calls
To find which bins have abnormally low/high coverage, `call.abnormal.cov` will derive P-value from the Z-score from a particular sample and use False Discovery Rate control. A simple bin merging strategy is natively available. In summary, looking at the Z-score of two consecutive bins, compared to the Z-scores of two randomly chosen bins, bins are joined or not. Other segmentation methods can be used on the Z-scores through `TODO` function.


### Running on computing clusters
`PopSV` can be used on a cluster using package `BatchJobs`. An example of an analysis using `BatchJobs` can
be found in folder `scripts`.

`BatchJobs` is a potent package but basic functions are enough in our situation. Here is a quick practical summary of `BatchJobs` commands used in the script:
* `makeRegistry` creates a registry used to manipulate jobs for a particular analysis step.
* `batchMap` adds jobs to a registry. Simply, the user gives a function and a list of parameters. One job per parameter will be created to compute the output of the function using this specific parameter.
* `submitJobs` submits the jobs to the cluster. This is where the queue, ,maximum computation time, number of core can be specified. Moreover, if needed, a subset of the jobs can be sent to the cluster. Functions `findNotDone` and `findErrors` are particularly useful to find which the jobs that didn't finish or were lost in the limbo of the cluster management process.
* `showStatus` outputs the status of the computations.
* `loadResult` retrieves the output of one specific job, while `reduceResultsList` retrieves output for all jobs into a list format.

Another important point about `BatchJobs` is its configuration for the computing cluster. An example of the configuration files can be found in the `scripts` folder:
* If present in the working directory, `.BatchJobs.R` is loaded when the `BatchJobs` package is loaded. It defines which template to use and `BatchJobs` functions. In practice, it loads another R script file (here `makeClusterFunctionsAdaptive.R`) with the functions to use. In `.BatchJobs.R` users would only need to change the email address where to send the log messages to.
* In `makeClusterFunctionsAdaptive.R`, users just need to check/replace `qsub`/`qdel`/`qstat` calls with the correct bash commands (sometimes `msub`/`canceljob`/`showq`). This file should also be in the working directory when `BatchJobs` is loaded.
* Finally `cluster.tmpl` is a template form of a job bash script that would be send to the cluster. There the correct syntax for the resources or parameters of the cluster are defined. This file should also be in the working directory when `BatchJobs` is loaded.
