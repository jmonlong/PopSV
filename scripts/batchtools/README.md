## Running PopSV in a HPC using *batchtools*

Place in the working directory: 

- `automatedPipeline-batchtools.R` which contains the pipeline functions.
- `batchtools.conf.R`, the configuration file for *batchtools*.
- `batchtools.slurm.tmpl`, the template for a job in your HPC.

Then in a R script, call the functions like in `run-PopSV-batchtools-automatedPipeline.R` (or copy and edit this script).

## Configuring your HPC

### Slurm

Slurm is used by Compute Canada on **Cedar**.

`batchtools.conf.R` is currently set up for Slurm (it calls *makeClusterFunctionsSlurm*).
When running the pipeline functions, this configuration file will be loaded and look for a template named `batchtools.slurm.tmpl`.

To work on Slurm (or at least on Cedar) you will have to specify your "account".
As you can see, there is a field for "account" in the template.
You can either replace the `<%= resources$account %>` by your account ID in the template as you would do when sending jobs manually, or specify your account when calling the R functions later.
For example, in `run-PopSV-batchtools-automatedPipeline.R`, we used `other.resources=list(account='rrg-bourqueg-ad')` to specify which account we want Cedar to use when sending jobs.

### TORQUE

TORQUE is used at McGill's Genome Center on **Abacus**.

To use TORQUE, `batchtools.conf.torque.R` should be renamed to `batchtools.conf.R` in the working directory.
This configuration will look for the `batchtools.torque.tmpl` template file when the pipeline is run.

### Others

To use other HPC systems you will need different configuration and template files.
There are a few examples of templates on [the *batchtools* GitHub](https://github.com/mllg/batchtools/tree/master/inst/templates).
In the configuration it about which *makeClusterFunctions* to use (more details in [*batchtools* documentation](https://mllg.github.io/batchtools/articles/batchtools.html)):  *makeClusterFunctionsDocker*, *makeClusterFunctionsInteractive*, *makeClusterFunctionsLSF*, *makeClusterFunctionsMulticore*, *makeClusterFunctionsOpenLava*, *makeClusterFunctionsSGE*, *makeClusterFunctionsSSH*, *makeClusterFunctionsSlurm*, *makeClusterFunctionsSocket*, *makeClusterFunctionsTORQUE*.

### Testing with a simple job

It's a good idea to check that everything is configured properly before trying to run the pipeline.
To test that sending jobs works you could try running the following commands:

```r
library(batchtools)

## To start again from scratch, manually remove the 'test' folder.
reg <- makeRegistry('test', seed=123)
## reg = loadRegistry('test', writeable=TRUE) ## If the registry has already been created before

test.f <- function(ii){
	return(mean(rnorm(10,ii)))
}

batchMap(reg=reg, test.f, 1:2)
submitJobs(reg=reg, ids=findJobs(reg=reg), resources=list(walltime='10:00', cores=1))
waitForJobs(reg=reg, sleep=10)
getStatus(reg=reg)
reduceResultsList(reg=reg)
```

These commands load the package, create a registry called *test*, define a function that will be run in the job, setup two jobs with this function and inputs 1 and 2, submit the jobs with a 10min walltime and 1 core per job, wait for the jobs to finish, show a summary of the status, list the results.
