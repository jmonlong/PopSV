***batchtools*  should be preferred over *BatchJobs*. See scripts in `scripts/batchtools`.** 

## Running PopSV in a HPC using *BatchJobs*

Place in the working directory: 

- `automatedPipeline-BatchJobs.R` which contains the pipeline functions.
- `BatchJobs_profile.R`, the configuration file for *BatchJobs*.
- `abacus.tmpl`, a template for a job in your HPC (here Abacus at McGill's Genome Center).
- `makeClusterFunctionsAdaptive`, the function to read the template file.

Then in a R script, call the functions like in `run-PopSV-batchjobs-automatedPipeline.R` (or copy and edit this script).

## Configuring your HPC

We used *BatchJobs* on TORQUE HPC only. 
Other HPC systems could be configured (see [BatchJobs doc](https://github.com/tudo-r/BatchJobs)) but we recommend using *batchtools* instead with the configuration in `scripts/batchjobs`.

TORQUE is used at McGill's Genome Center on **Abacus** and by Compute Canada on **Guillimin**.

To configure *BatchJobs* from these files, you should:

1. Check the name of the template file called in `BatchJobs_profile.R`. Currently `guillimin.tmpl`. Change to `abacus.tmpl` for Abacus.
1. Check the commands in `makeClusterFunctionsAdaptive.R`. Currently `qsub`/`canceljob`/`qstat`. Some HPC prefer `msub`, `showq`, `qdel`, etc.

When running the pipeline on Guillimin, add `other.resources=list(supervisor.group="XXXXX")` with the account to use. 
This should be added as an argument of the *autoGCcounts*, *autoExtra* and *autoNormTest* functions.
