---
layout: page
title: Cluster management in R
---

## *BatchJobs* package

*BatchJobs* is a potent package to communicate with a cluster, i.e. send jobs, check their status, eventually rerun, retrieve the results.
PopSV has been designed into separate steps to be ran on a computing cluster using *BatchJobs*.

A one-command version is also available there. It's basically a wrapper for the basic analysis steps with some useful functions (running custom steps, stop/restart). However it's less flexible.

## Installation and configuration

The package can be installed through CRAN

```r
install.packages("BatchJobs")
```

The most important part about *BatchJobs* is its configuration for your computing cluster. It's not long but should be done carefully. And once this is done correctly, the rest follows nicely.

You will need to create **3 files** : a cluster template, parser functions, and a `.BatchJobs.R` configuration file. I would recommend to put these 3 files **in the root of your personal space**, i.e. `~/`. You could put them in the project folder, but then it means you have to copy them each time you create/run another project. Putting them in your root means that they will be used by default by *BatchJobs*.

### Cluster template

A cluster template is a template form of a job bash script that you would send through `qsub`/`msub`. There you define the placeholder for the resources or parameters of the job. This file will be parsed by *BatchJobs*. 

For example, we create a `guillimin.tmpl` file for our cluster [Guillimin](http://www.hpc.mcgill.ca/) like this :

```sh
#PBS -N <%= job.name %>
#PBS -j oe
#PBS -o <%= log.file %>
#PBS -l walltime=<%= resources$walltime %>
#PBS -l nodes=<%= resources$nodes %>:ppn=<%= resources$cores %>
#PBS -A <%= resources$supervisor.group %>
#PBS -q <%= resources$queue %>
#PBS -V

## Run R:
## we merge R output with stdout from PBS, which gets then logged via -o option
R CMD BATCH --no-save --no-restore "<%= rscript %>" /dev/stdout
```

Placeholders are in the form of `<%= resources$walltime %>`. *BatchJobs* will insert there the element define by `walltime` element in the `resources` list (see later). Although you most likely won't have to change these placeholders, you might need to update the lines if your cluster uses a different syntax. 

### Parser functions

Parser functions are saved in a R script, called for example `makeClusterFunctionsAdaptive.R`. This will parse the template and create the actual commands to send, cancel and check jobs.

Likely you just need to check/replace `qsub`/`qdel`/`qstat` calls with the correct bash commands (sometimes `msub`/`canceljob`/`showq`).

From our file, these are the lines you might need to change :

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

In practice, **you won't have to write this part** as we provide full pipelines. You might still need to change a bit the resources of the jobs (they might change from one cluster to another). More precisely I'm talking about the `resource=` parameter in the `submitJobs` command. After doing this, if you are not interested in more details, you can jump directly to the next section for an overview of a pipeline script.

Here is a quick summary of *BatchJobs* commands used in the scripts:

* `makeRegistry` creates a registry used to manipulate jobs for a particular analysis step.
* `batchMap` adds jobs to a registry. Simply, you give it a function and a list of parameters. One job per parameter will be created to compute the output of the function using this specific parameter.
* `submitJobs` submits the jobs to the cluster. This is where the queue, maximum computation time, number of cores can be specified. Moreover, if needed, a subset of the jobs can be sent to the cluster. Functions `findNotDone` and `findErrors` are particularly useful to find which jobs didn't finish or were lost in the limbo of the cluster management process.
* `showStatus` outputs the status of the computations.
* `loadResult` retrieves the output of one specific job, while `reduceResultsList` retrieves output for all jobs into a list format.

```r
getBC.reg <- makeRegistry(id="getBC")
getBC.f <- function(file.i, gc.reg, files.df){
  library(PopSV)
  bins.df = loadResult(gc.reg,1)
  bin.bam(files.df$bam[file.i], bins.df, files.df$bc[file.i])
}
batchMap(getBC.reg, getBC.f,1:nrow(files.df), more.args=list(gc.reg=getGC.reg, files.df=files.df))
submitJobs(getBC.reg, findNotDone(getBC.reg), resources=list(walltime="20:0:0", nodes="1", cores="1"), wait=function(retries) 100, max.retries=10)
```

Here we want to get the bin counts of each sample. We create a registry called *getBC*. Then the function that will get the bin counts for a sample. Here `file.i` is the index of the sample and will be different for each job sent by *BatchJobs*. The other parameters are common to all jobs. Within the function, we load the package and useful data and run the function we want. `batchMap` creates a job per sample id and link the function we just defined. The jobs are finally submitted to the cluster with specific number of cores, wall time, ... 

We can check the status of the job with `showStatus(getBC.reg)` command. 

## Pipeline worflow

The general idea is to have one script per analysis (e.g. bin size, project). Examples of pipeline scripts can be found in the `scripts` folder of the GitHub repository.

The script doesn't need to be ran each time from the start but rather ran step by step, most likely at separate times. Think about R as a new shell: in R the status of the jobs in the clusters are checked, rerun, etc. Nothing will be directly computed by this script so it can be ran in a login node.

When one step sends jobs to the cluster through *BatchJobs*, the user can exit R, log out, have a coffee, think about all the time saved thanks to *BatchJobs* and then open everything again and continue. No need to rerun everything, just load the libraries and the registry of the steps to check and continue.
