# Cluster configuration with BatchJobs

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
