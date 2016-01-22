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

#### Count 
male.GCcounts = autoGCcounts("filesMale.RData", "bins.RData", file.suffix="male", lib.loc="/home/jmonlong/R/PopSVforPaper", status=TRUE)
female.GCcounts = autoGCcounts("filesFemale.RData", "bins.RData", file.suffix="female", lib.loc="/home/jmonlong/R/PopSVforPaper", skip=1, status=TRUE)

#### Run PopSV
resFemale.df = autoNormTest("filesFemale.RData", "bins.RData", file.suffix="female", lib.loc="/home/jmonlong/R/PopSVforPaper")
resMale.df = autoNormTest("filesMale.RData", "bins.RData", file.suffix="male", lib.loc="/home/jmonlong/R/PopSVforPaper")

#### Merge results
res.df = rbind(resFemale.df, resMale.df)
res.filt.df = sv.summary.interactive(res.df) ## Run locally because it opens an interactive web browser application
