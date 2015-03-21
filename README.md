PopSV
=====
PopSV is a structural variation (SV) detection method from high-throughput sequencing. 
Abnormal Read-Depth signal is detected by using a population of samples ase reference. Thanks to this population
view the whole genome can be robustly interrogated, including regions of low mappability. Moreover, any divergence from
the reference samples are detected, even if the signal is incomplete, e.g. tumoral aberrations or SV involving repeats.

**Warning: PopSV package is still in active development.**

### Installation
To install the latest development version: `devtools::install_github("jmonlong/PopSV")`. This requires `devtools` package (more information [here](https://github.com/hadley/devtools)) which can be easily installed with `install.packages("devtools")`. 

Some [Bioconductor](http://bioconductor.org/) packages are also necessary and not installed automatically. Running the following should be sufficient :
```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19", "Rsamtools")
```

**R 3.1 or higher** is required.


### Analysis steps
Examples of a analysis, for local computation or computing cluster usage, can be found on the `scripts` folder. For more information on a specific function, see the manual or access the documentation through `?the.function.name`.

#### Input files
The analysis can start directly from the BAM files. Each BAM file needs to be **sorted** and **indexed** (see [samtools](http://www.htslib.org/)).

A tabular separated values (*tsv*) file with the name of the sample (in column named *sample*) and the path to the corresponding BAM (in column named *bam*) is imported and  given to the `init.filenames`. This function will create the path and file names for the different files created and used throughout the analysis.

Then, the regions of interest or *bins* have to be defined. `fragment.genome.hp19` automatically fragments hg19 genome into non-overlapping consecutive windows of a specified size. However PopSV can perform with any type of windows. It is still recommended to define non-overlapping windows and, for computation reasons, no more than a total of several million of bins. If a custom definition is used it should follow the BED format with columns named *chr*, *start* and *end*.

Finally, the GC content of each bin can be computed, for hg19, using function `getGC.hg19`. If another genome is to be used, GC content should be define following BedGraph format, i.e with columns named *chr*, *start*, *end* and *GCcontent*.

#### Counting reads
Reads are counted in each bin to measure coverage. `bin.bam` function will count reads in each bin for a given sample.

Eventually, this can be done externally, e.g. using [bedtools coverage](http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html). The final count file should have these four columns: *chr*, *start*, *end* and *bc* (for bin count).

#### GC bias correction
GC bias is corrected using a LOESS model. Using this model, a normalization coefficient is computed for each bin based on its GC content. This step is performed by `correct.GC` function.

#### Sample Quality Control
The last "pre-processing" step aims at defining the set of reference samples. These samples will define "normal" coverage. A natural set of reference samples are the controls in case/control studies or normal samples in normal/tumor paired samples designs. The more the better but these samples should also be homogeneous to get optimal detection power.

`qc.samples.summary` opens an interactive web-browser application to explore how homogeneous the samples are. The samples are represented and clustered using the first two principal components. The user decides how many clear clusters are visible and which one will be analyzed. **This step is usually not necessary** because samples analyzed should have been sequenced with similar protocol and technologies. It mostly a safety check. Here a subset of the bins can be used, e.g. using `quick.count` function.

`qc.samples` function will join all bin count files and produce some graphs on the set of reference samples. If too many reference samples are available (lucky you), `nb.ref.samples=` parameters can be used. 200 reference samples would be enough.

#### Normalization of the reference samples

This step is extremely important to avoid sample-specific bias being picked up as abnormal coverage. Several normalization approaches are available but, for now, targeted normalization (`tn.norm` function) should be used. Fortunately, it's the one that works best (in my experience). Other function are/will be `pca.norm` for Principal Component removal normalization, `quant.norm` for quantile normalization, `medvar.norm` for a more naive normalization of the coverage median and variance.

The output of these functions is a data.frame with the normalized bin counts as well as, sometimes, some metrics on the normalization efficiency.

#### Z-score computation for the reference samples

From the normalized bin counts the mean and standard deviation across reference samples is computed for each bin. A Z-score and fold-change estimate are derived for each reference sample. Function `z.comp` performs this step.

#### Testing other samples

At this point the rest of the samples can be tested using `tn.test.sample` function. Both normalization and Z-score computation are performed by the function.

#### Abnormal coverage calls
To find which bins have abnormally low/high coverage, `call.abnormal.cov` will derive P-value from the Z-score from a particular sample and use False Discovery Rate control. A stitching bin merging strategy is available and recommended through parameter `merge.cons.bins="stitch"`. Other parameters are described in the function documentation, such as the strategy for P-value definition or extra tricks (usually for cancer samples, e.g. aneuploid chromosome removal). 

#### Visualizing the results
Function `sv.summary.interactive` opens an interactive web-browser application to look at the results. The number of calls across samples, distribution of copy-number estimates or frequency is visualized. The user can play with different stringency filters in order to retrieve the ones that gives the best result quality.


### Running on computing clusters
`PopSV` can be used on a cluster using package `BatchJobs`. An example of an analysis using `BatchJobs` can be found in folder `scripts`.

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

The general idea is to have one script, such as the one in the `scripts` folder, per analysis (e.g. bin size, project). The script should not be run each time from the start but rather ran step by step, likely at separate times. Think as R as the new shell: in R the status of the jobs in the clusters are checked, rerun, etc. Indeed when one step sends jobs to the cluster through `BatchJobs`, the user can exit R, logout, have a coffee, think about all the time saved thanks to `BatchJobs` and then open everything again and continue. No need to rerun everything, just load the libraries and the registry of the steps to check and continue.
