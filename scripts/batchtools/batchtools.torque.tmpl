#!/bin/bash
<%
## Check some resources and set sane defaults
resources$nodes = if (is.null(resources$nodes)) 1L else resources$nodes
-%>

#PBS -N <%= job.name %>
#PBS -o <%= log.file %>
#PBS -l walltime=<%= resources$walltime %>
#PBS -l nodes=<%= resources$nodes %>:ppn=<%= resources$cores %>
#PBS -l vmem=<%= resources$memory %>
#PBS -j oe
<%= if (array.jobs) sprintf("#PBS -t 1-%i", nrow(jobs)) else "" %>

## export value of DEBUGME environemnt var to slave
export DEBUGME=<%= Sys.getenv("DEBUGME") %>

## Run R
Rscript -e 'batchtools::doJobCollection("<%= uri %>")'
