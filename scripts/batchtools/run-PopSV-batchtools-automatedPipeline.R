library(PopSV)
source('automatedPipeline-batchtools.R')

bam.files = read.table("bams.tsv", as.is=TRUE, header=TRUE)
bin.size = 1e3

#### Init file names and construct bins
files.df = init.filenames(bam.files, code="example")
bins.df = fragment.genome.hg19(bin.size)
save(files.df, file="files.RData")
save(bins.df, file="bins.RData")
####

res.GCcounts = autoGCcounts("files.RData", "bins.RData", other.resources=list(account='rrg-bourqueg-ad'))

res.forQC = autoExtra("files.RData", "bins.RData", do=1, other.resources=list(account='rrg-bourqueg-ad')))
qc.samples.cluster(res.forQC) ## Run locally because it opens an interactive web browser application

res.df = autoNormTest("files.RData", "bins.RData", other.resources=list(account='rrg-bourqueg-ad')))

res.filt.df = sv.summary.interactive(res.df) ## Run locally because it opens an interactive web browser application

