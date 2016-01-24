#### Installation
## devtools::install_github("jmonlong/PopSV")

library(BatchJobs)
library(PopSV)
source("automatedPipeline.R")

bam.files = read.table("bam-samples.tsv", as.is=TRUE, header=TRUE)
bin.size = 5e3

#### Init file names and construct bins
bins.df = fragment.genome.hg19(bin.size, XY.chr=TRUE)
save(bins.df, file="bins.RData")
files.df = init.filenames(subset(bam.files, gender=="male"), code="5kbpMale")
save(files.df, file="filesMale.RData")
files.df = init.filenames(subset(bam.files, gender=="female"), code="5kbpFemale")
save(files.df, file="filesFemale.RData")
####

## Reference samples: normal samples
ref.samples = subset(bam.files, status=="normal")$sample

## Count
male.GCcounts = autoGCcounts("filesMale.RData", "bins.RData", file.suffix="male")
female.GCcounts = autoGCcounts("filesFemale.RData", "bins.RData", file.suffix="female", skip=1)
##


#### Run PopSV
resMale.df = autoNormTest("filesMale.RData", "bins.RData", file.suffix="male", ref.samples=ref.samples)
resFemale.df = autoNormTest("filesFemale.RData", "bins.RData", file.suffix="female", ref.samples=ref.samples)


#### Merge results
res.df = rbind(resFemale.df, resMale.df)
## Eventually merge useful information on the samples
res.df = merge(res.df, bam.files[,c("sample","gender","status")])
res.filt.df = sv.summary.interactive(res.df) ## Run locally because it opens an interactive web browser application

