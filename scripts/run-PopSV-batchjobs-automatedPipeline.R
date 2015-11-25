#### Installation
## devtools::install_github("jmonlong/PopSV")

library(BatchJobs)
library(PopSV)

bam.files = read.table("bams.tsv", as.is=TRUE, header=TRUE)
bin.size = 1e3

#### Init file names and construct bins
files.df = init.filenames(bam.files, code="example")
bins.df = fragment.genome.hp19(bin.size)
save(files.df, file="files.RData")
save(bins.df, file="bins.RData")
####

source("automatedPipeline.R")

res.GCcounts = autoGCcounts("files.RData", "bins.RData")

qc.samples.cluster(res.GCcounts) ## Run locally because it opens an interactive web browser application

res.df = autoNormTest("files.RData", "bins.RData")

sv.summary.interactive(res.df) ## Run locally because it opens an interactive web browser application

