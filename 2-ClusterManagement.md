---
layout: pagetoc
title: Cluster management in R
permalink: /2-ClusterManagement.md/
---

*batchtools* is the successor of *BatchJobs* and will be used from now on as it is more stable and more flexible. 
To see the old version of this page, detailing how *BatchJobs* was used to run PopSV, see [this page](https://github.com/jmonlong/PopSV/blob/master/2-ClusterManagement-BatchJobs.md).

## *batchtools* package

[*batchtools*](https://mllg.github.io/batchtools/) is a useful package to communicate with a computing cluster: send jobs, check their status, eventually rerun them, retrieve the results.
PopSV has been designed into separate steps to run more easily on a computing cluster using *batchtools*. 
Thanks to the multi-step workflow, the computation is parallelized as much as possible (sometimes by sample, other times by genomic regions).

Instead of running each step manually, we recommend using the two-commands wrappers (see *Automated run* below). 
It's basically a wrapper for the basic analysis steps with some useful functions (running custom steps, stop/restart). 
It should be sufficient for most analysis but it's less flexible and if you want to change some specific parameters you might have to tweak the code within the functions.

## Installation and configuration

*batchtools* package can be installed through [CRAN](https://www.cran.r-project.org/):

```r
install.packages("batchtools")
```

The most important step is **configuring it for your computing cluster**. 
It's not long to do and once this is done correctly, the rest follows nicely.

You will need to place **2 files** in the working directory of your project:

+ a cluster template
+ a configuration file.

I put some **examples of configuration and template files for Slurm and TORQUE HPC** in the [`scripts` folder of the GitHub repo](https://github.com/jmonlong/PopSV/tree/master/scripts/batchtools).
If you are using Slurm/TORQUE or some of the HPC from Compute Canada (e.g. Cedar/Guillimin), you can use these files with only minimal edits (see instructions on GitHub).

### Cluster template

A cluster template is a template form of the bash script that you would send through `qsub`/`msub`. 
There you define the placeholder for the resources or parameters of the job. This file will be parsed by *batchtools*.

For the Compute Canada cluster [Cedar](https://docs.computecanada.ca/wiki/Cedar), I created use a template that looks like this:

```sh
#!/bin/bash

#SBATCH --job-name=<%= job.name %>
#SBATCH --output=<%= log.file %>
#SBATCH --error=<%= log.file %>
#SBATCH --time=<%= resources$walltime %>
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<%= resources$cores %>
#SBATCH --mem-per-cpu=<%= resources$memory %>
#SBATCH --account=<%= resources$account %>
<%= if (!is.null(resources$partition)) sprintf(paste0("#SBATCH --partition='", resources$partition, "'")) %>
<%= if (array.jobs) sprintf("#SBATCH --array=1-%i", nrow(jobs)) else "" %>

## Initialize work environment like
## source /etc/profile
## module add ...

## Export value of DEBUGME environemnt var to slave
export DEBUGME=<%= Sys.getenv("DEBUGME") %>

## Run R:
## we merge R output with stdout from SLURM, which gets then logged via --output option
Rscript -e 'batchtools::doJobCollection("<%= uri %>")'
```

Placeholders are in the form of `<%= resources$walltime %>`. 
*batchtools* will insert there the value defined by `walltime` element in the `resources` list (see later in `submitJobs` command). 
Although you most likely won't have to change these placeholders, you might need to update the lines if your cluster uses a different syntax. 
For example, in our cluster, we need to give an account ID for our lab, that will be passed in the template by `<%= resources$account %>`.

In order to easily use the pipelines provided with PopSV package, you should **use at least the placeholders `walltime`, `cores` and `nodes`**.
For the other arguments (e.g. queue, lab ID) you can either hard-code them in the template or pass it in the R pipeline functions with `other.resources=list(...)`.

### Configuration file

The `batchtools.conf.R` file will be read when running the R pipeline functions.
It will load the template file defined above. 
There is nothing special to change here. 
Just making sure that the command fits the HPC:

- For Slurm, `batchtools.conf.R` contains `cluster.functions = makeClusterFunctionsSlurm()` and will look for a template file called `batchtools.slurm.tmpl`.
- For TORQUE, `batchtools.conf.R` contains `cluster.functions = makeClusterFunctionsTORQUE()` and will look for a template file called `batchtools.torque.tmpl`.

## Running the pipeline

### Automated run

Two wrapper functions around *batchtools* allows you to run PopSV without manually sending the jobs for each steps. 
These two functions (`autoGCcounts` and `autoNormTest`) are located in [`automatedPipeline-batchtools.R`](https://github.com/jmonlong/PopSV/tree/master/scripts/batchtools). 
A full analysis can be run like this:

```r
## Load package and wrapper
library(PopSV)
source("automatedPipeline-batchtools.R")
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

The advantage of this wrapper is a easier management of the cluster and pipeline. 
However it's not so flexible: if a step need to be changed for some reason, you might have to change it within the `automatedPipeline-batchtools.R` script.

Still, a few parameters can be passed to the two functions for the user's convenience:

+ Use `lib.loc=` if you installed PopSV in a specific location. The value will be passed to `library(PopSV)`.
+ `redo=` can be used to force a step to be redone (i.e. previous jobs deleted and re-submitted). E.g. `redo=5` to redo step 5.
+ `other.resources=` to specify resources for the jobs to match the template (see HPC configuration section above). We use this to specify queues/accounts when the HPC requires it.
+ `resetError=TRUE` to reset jobs that had errors and rerun them. Better than a *redo* because the jobs that are done don't need to be rerun.
+ `rewrite=TRUE` will force the normalized bin counts and normalization stats to be rewritten.
+ `file.suffix=` to add a suffix to the temporary files. This is useful when the pipeline is run several times on the same folder, for example when splitting the samples in batches (e.g. presence of batch effects, male/female split for XY chrs).
+ `step.walltime=` the walltime for each step. See in the `automatedPipeline-batchtools.R` script for default values. 
+ `step.cores=` the number of cores for each step. See in the `automatedPipeline-batchtools.R` script for default values. 
+ `status=TRUE` will print the status of the jobs for each steps and the end of the log of jobs with errors.
+ `skip=` to skip some steps.


As an example, the `run-PopSV-batchtools-automatedPipeline.R` script in [the scripts folder](https://github.com/jmonlong/PopSV/tree/master/scripts/batchtools) shows how PopSV is run using these wrappers on Cedar (Slurm HPC). 

### Practical details

- The configuration files and `automatedPipeline-batchtools.R` script should be in the working directory. 
- Use different `file.suffix` if PopSV is run several times in the same folder (e.g. different bin size, sample batches).
- The master script (~`run-PopSV-batchtools-automatedPipeline.R`) can be left running on a login node of the HPC because it doesn't compute anything, it just sends jobs and wait. Even better, let it run on a [*screen*](https://www.gnu.org/software/screen/manual/screen.html) so that you can detach it and disconnect from the server without stopping the pipeline.
- The paths and folder structure is saved in the `files.df` data.frame, originally created by  `init.filenames` function. 

## Preparing and submitting a job with *batchtools* 

In practice, **you don't have to write this part**, it's what the `automatedPipeline-batchtools.R` functions are made of.
If ever you want to tweak these functions or just use *batchtools* for something else, here is how we use it.

Useful functions:

- `loadRegistry` creates a registry used to manipulate jobs for a particular analysis step. Use `writable=TRUE` if the registry already exists. 
- `batchMap` adds jobs to a registry. You give it a function and a list of parameters. One job per parameter will be created to compute the output of the function using this specific parameter. `more.args=` to provide additional arguments (same for all jobs).
- `submitJobs` submits the jobs to the cluster. This is where the walltime time, number of cores, etc can be specified. Moreover, if needed, a subset of the jobs can be sent to the cluster. Functions `findNotDone` and `findErrors` are particularly useful to find which jobs didn't finish or were lost in the limbo of the cluster management process.
- `getStatus` outputs the status of the computations.
- `loadResult` retrieves the output of one specific job, while `reduceResultsList` retrieves output for all jobs into a list format.
- `waitForJobs` waits for all the jobs to finish.

For example, to run the step that retrieve the bin count in each BAM files, a simple version would look like this:

```r
reg <- loadRegistry("getBC", writeable=TRUE)
getBC.f <- function(file.i, bins.f, files.df){
  library(PopSV)
  load(bins.f)
  bin.bam(files.df$bam[file.i], bins.df, files.df$bc[file.i])
}
batchMap(reg=reg, getBC.f,1:nrow(files.df), more.args=list(bins.f="bins.RData", files.df=files.df))
submitJobs(reg=reg, findNotDone(reg=reg), resources=list(walltime="20:0:0", nodes="1", cores="1"))
getStatus(reg=reg)
```

Here we want to get the bin counts of each sample. 
We create a registry called *getBC*. 
Then we define the function that will get the bin counts for a sample. 
The first parameter of this function (here `file.i` which is the index of the sample) will be different for each job sent by *batchtools*.
Other parameters are the same in all jobs. 
Within the function, we load the package and useful data and run the instructions we want. 
`batchMap` creates one job per sample ID and links the function we've just defined. 
The jobs are finally submitted to the cluster with the desired number of cores, walltime, etc
we can check the status of the job with `getStatus` function.
