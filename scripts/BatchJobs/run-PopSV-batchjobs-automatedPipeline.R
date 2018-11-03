#### Installation
## devtools::install_github("jmonlong/PopSV")

library(BatchJobs)
loadConfig('BatchJobs_profile.R')

library(PopSV)

bam.files = read.table("bams.tsv", as.is=TRUE, header=TRUE)
bin.size = 1e3

#### Init file names and construct bins
files.df = init.filenames(bam.files, code="example")
bins.df = fragment.genome.hg19(bin.size)
save(files.df, file="files.RData")
save(bins.df, file="bins.RData")
####

source("automatedPipeline-BatchJobs.R")

res.GCcounts = autoGCcounts("files.RData", "bins.RData")

res.forQC = autoExtra("files.RData", "bins.RData", do=1)
qc.samples.cluster(res.forQC) ## Run locally because it opens an interactive web browser application

res.df = autoNormTest("files.RData", "bins.RData")

res.filt.df = sv.summary.interactive(res.df) ## Run locally because it opens an interactive web browser application

