---
layout: page
title: Cluster management in R
permalink: /2-ClusterManagement.md/
---

## *BatchJobs* package

[*BatchJobs*](https://github.com/tudo-r/BatchJobs) is a potent package to communicate with a cluster, i.e. send jobs, check their status, eventually rerun, retrieve the results.
PopSV has been designed into separate steps to be ran more easily on a computing cluster using *BatchJobs*.

A two-commands version is also available (see *Automated run* below. It's basically a wrapper for the basic analysis steps with some useful functions (running custom steps, stop/restart). However it's less flexible.

## Installation and configuration

*BatchJobs* package can be installed through [CRAN](https://www.cran.r-project.org/)

```r
install.packages("BatchJobs")
```

The most important step is **configuring it for your computing cluster**. It's not long but should be done carefully. And once this is done correctly, the rest follows nicely.

You will need to create **3 files** :

+ a cluster template
+ a script with parser functions
+ a `.BatchJobs.R` configuration file.

I would recommend to put these 3 files **in the root of your personal space**, i.e. `~/`. You could put them in the project folder, but then it means you have to copy them each time you create/run another project. Putting them in your root means that they will always be used by default by *BatchJobs*.

### Cluster template

A cluster template is a template form of the bash script that you would send through `qsub`/`msub`. There you define the placeholder for the resources or parameters of the job. This file will be parsed by *BatchJobs*.

For our cluster [Guillimin](http://www.hpc.mcgill.ca/), we create a `guillimin.tmpl` file like this :

```sh
#PBS -N <%= job.name %>
#PBS -j oe
#PBS -o <%= log.file %>
#PBS -l walltime=<%= resources$walltime %>
#PBS -l nodes=<%= resources$nodes %>:ppn=<%= resources$cores %>
#PBS -A bws-221-ae
#PBS -V

## Run R:
## we merge R output with stdout from PBS, which gets then logged via -o option
R CMD BATCH --no-save --no-restore "<%= rscript %>" /dev/stdout
```

Placeholders are in the form of `<%= resources$walltime %>`. *BatchJobs* will insert there the value defined by `walltime` element in the `resources` list (see later in `submitJobs` command). Although you most likely won't have to change these placeholders, you might need to update the lines if your cluster uses a different syntax. For example, in our cluster, we need to give an ID for our lab with `-A`.

In order to easily use the pipelines provided with PopSV package, I would recommend to **put exactly the placeholders** `walltime`, `cores` and `nodes` but to hard-code the rest (e.g. queue, lab ID) in the template.

### Parser functions

Parser functions are saved in a R script, called for example `makeClusterFunctionsAdaptive.R`. This will parse the template and create the actual commands to send, cancel and check jobs.

Likely you just need to check/replace `qsub`/`qdel`/`qstat` calls with the correct bash commands (sometimes `msub`/`canceljob`/`showq`).

From [our file](https://github.com/jmonlong/PopSV/blob/forPaper/scripts/makeClusterFunctionsAdaptive.R), these are the lines you might need to change :

```sh
res =  BatchJobs:::runOSCommandLinux("qsub", outfile, stop.on.exit.code = FALSE)
cfKillBatchJob("canceljob", batch.job.id)
BatchJobs:::runOSCommandLinux("showq", "-u $USER")$output
```

### `.BatchJobs.R` configuration file
`.BatchJobs.R`  is just the configuration files that links the two other files. You don't really need to change it. Eventually change the email address.

It looks like this :

```r
source("~/makeClusterFunctionsAdaptive.R")
cluster.functions <- makeClusterFunctionsAdaptive("~/guillimin.tmpl")
mail.start <- "none"
mail.done <- "none"
mail.error <- "none"
mail.from <- "<jean.monlong@mail.mcgill.ca>"
mail.to <- "<jean.monlong@mail.mcgill.ca>"
```

*Note: If `.BatchJobs.R` files are present at both `~/` and the project folder, the one in the project folder will override the parameters.*


## Sending Jobs

In practice, **you won't have to write this part** as we provide full pipelines. You might still need to **change a bit the resources of the jobs** (they might change from one cluster to another). More precisely I'm talking about the `resource=` parameter in the `submitJobs` command. After doing this, if you are not interested in more details, you can jump directly to the next section for an overview of a pipeline script.

Otherwise here is a quick summary of *BatchJobs* commands used in the scripts:

* `makeRegistry` creates a registry used to manipulate jobs for a particular analysis step.
* `batchMap` adds jobs to a registry. Simply, you give it a function and a list of parameters. One job per parameter will be created to compute the output of the function using this specific parameter.
* `submitJobs` submits the jobs to the cluster. This is where the queue, maximum computation time, number of cores can be specified. Moreover, if needed, a subset of the jobs can be sent to the cluster. Functions `findNotDone` and `findErrors` are particularly useful to find which jobs didn't finish or were lost in the limbo of the cluster management process.
* `showStatus` outputs the status of the computations.
* `loadResult` retrieves the output of one specific job, while `reduceResultsList` retrieves output for all jobs into a list format.

For example, to run the step that retrieve the bin count in each BAM files it looks like this :

```r
getBC.reg <- makeRegistry(id="getBC")
getBC.f <- function(file.i, bins.f, files.df){
  library(PopSV)
  load(bins.f)
  bin.bam(files.df$bam[file.i], bins.df, files.df$bc[file.i])
}
batchMap(getBC.reg, getBC.f,1:nrow(files.df), more.args=list(bins.f="bins.RData", files.df=files.df))
submitJobs(getBC.reg, findNotDone(getBC.reg), resources=list(walltime="20:0:0", nodes="1", cores="1"))
showStatus(getBC.reg)
```

Here we want to get the bin counts of each sample. We create a registry called *getBC*. Then we define the function that will get the bin counts for a sample. The first parameter of this function (here `file.i` which is the index of the sample) will be different for each job sent by *BatchJobs*. Other parameters are common to all jobs. Within the function, we load the package and useful data and run the instructions we want. `batchMap` creates one job per sample ID and links the function we've just defined. The jobs are finally submitted to the cluster with the desired number of cores, wall time, etc

We can check the status of the job with `showStatus(getBC.reg)` command.

## Pipeline workflow

### Automated run

Two wrapper functions around *BatchJobs* allows you to run PopSV without manually sending the jobs for each steps. These two functions (`autoGCcounts` and `autoNormTest`) are located in [`automatedPipeline.R`](https://github.com/jmonlong/PopSV/blob/forPaper/scripts/automatedPipeline.R). Now, a full analysis can be run like this:

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

The advantage of this wrapper is a easier management of the cluster and pipeline. However it's not so flexible: if a step need to be changed for some reason, you might have to change it within the `automatedPipeline.R` script.

Still, a few parameters can be passed to the two functions for the user convenience:

+ Use `lib.loc=` if you installed PopSV in a specific location. The value will be passed to `library(PopSV)`.
+ `redo=` can be used to force a step to be redone (i.e. previous jobs deleted and re-submitted). E.g. `redo=5` to redo step 5.
+ `rewrite=TRUE` will force the normalized bin counts and normalization stats to be rewritten.
+ `file.suffix=` to add a suffix to the temporary files. This is useful when the pipeline is run several times on the same folder, for example when splitting the samples in batches (e.g. presence of batch effects, male/female split for XY chrs).

We provide a [script](https://github.com/jmonlong/PopSV/blob/forPaper/scripts/run-PopSV-batchjobs-automatedPipeline.R) to PopSV using these wrappers. In addition, when we want to analyze X and Y chromosomes, the samples have to be split and these wrappers come handy to run easily three analysis (see this [example](https://github.com/jmonlong/PopSV/blob/forPaper/scripts/run-PopSV-XY-batchjobs-automatedPipeline.R)).

### Step-by-step manual run

The general idea is to have one script per analysis (e.g. bin size, project). Each such analysis should be in its own folder to avoid possible confusion between temporary files. Examples of pipeline scripts can be found in the [`scripts` folder of the GitHub repository](https://github.com/jmonlong/PopSV/tree/forPaper/scripts).

Because it manipulates large data (BAM files, genome-wide coverage) and large sample sizes, PopSV was designed to create and work with intermediate files. The management of these files are mostly handle automatically. In practice all the important path and folder structure is saved in the `files.df` data.frame, originally created by  `init.filenames` function. For this reason, the results of each analysis steps are saved as the local files so that the next steps can be run later without the need for you to think about what to save etc.

So the script doesn't need to be ran each time from the start but rather ran step by step. In practice you often have to wait a couple of hours for some step to compute. Think about R as a new shell: you would open R, check the status of the jobs in the clusters, rerun them if necessary, or start the next step, etc. You can run R on this master script in a login node because nothing will be directly computed there, the real computation are sent as actual jobs.

After one step sends jobs to the cluster, the user can exit R, log out, have a coffee, think about all the time saved thanks to *BatchJobs* and then open everything again later and continue. No need to rerun everything, just load the libraries and the registries (e.g. running `getBC.reg <- makeRegistry(id="getBC")` again) you want to check.
