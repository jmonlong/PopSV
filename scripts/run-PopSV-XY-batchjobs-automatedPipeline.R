#### Installation
## devtools::install_github("jmonlong/PopSV")

library(BatchJobs)
library(PopSV)

bam.files = read.table("bams.tsv", as.is=TRUE, header=TRUE)
bin.size = 1e3

#### Init file names and construct bins
files.df = init.filenames(bam.files, code="example")
save(files.df, file="files.RData")
bins.df = fragment.genome.hg19(bin.size, XY.chr=TRUE)
save(bins.df, file="bins.RData")
####

#### Create and save input for the 3 runs
## All samples, autosomes
files.df = init.filenames(bam.df, code="example")
save(files.df, file="files.RData")
load("bins.RData")
bins.df = subset(bins.df, chr %in% 1:22)
save(bins.df, file="binsAll.RData")
## Females, chr X
files.df = init.filenames(subset(bam.df, ped=="Mother"), code="exampleFemale")
save(files.df, file="filesFemale.RData")
load("bins.RData")
bins.df = subset(bins.df, chr %in% "X")
save(bins.df, file="binsFemale.RData")
## Males
files.df = init.filenames(subset(bam.df, ped=="Father"), code="exampleMale")
save(files.df, file="filesMale.RData")
load("bins.RData")
bins.df = subset(bins.df, chr %in% c("X","Y"))
save(bins.df, file="binsMale.RData")
####

source("automatedPipeline.R")

#### Run PopSV
## All samples
res.GCcounts = autoGCcounts("files.RData", "binsAll.RData")
resAll.df = autoNormTest("files.RData", "binsAll.RData")
## Females
res.GCcounts = autoGCcounts("filesFemale.RData", "binsFemale.RData", file.suffix="female")
resFemale.df = autoNormTest("filesFemale.RData", "binsFemale.RData", file.suffix="female", lib.loc="/home/jmonlong/R/PopSVforPaper")
## Males
res.GCcounts = autoGCcounts("filesMale.RData", "binsMale.RData", file.suffix="male")
resMale.df = autoNormTest("filesMale.RData", "binsMale.RData", file.suffix="male")

#### Merge results
res.df = rbind(resAll.df, resFemale.df, resMale.df)
res.filt.df = sv.summary.interactive(res.df) ## Run locally because it opens an interactive web browser application
