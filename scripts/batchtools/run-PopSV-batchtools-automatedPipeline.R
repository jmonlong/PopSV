library(PopSV)
source('automatedPipeline-batchtools.R')

genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
## or
## genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

##
## Preparation
## Run only once to create the files files.RData and bins.RData
##

bam.files = read.table("bams.tsv", as.is=TRUE, header=TRUE)
bin.size = 1e3

#### Init file names and construct bins
files.df = init.filenames(bam.files, code="example")
bins.df = fragment.genome(bin.size, genome=genome)
save(files.df, file="files.RData")
save(bins.df, file="bins.RData")
####


##
## Analysis
## Can be stopped and restarted. No need to rerun the preparation commands
##

## Bin and count reads in each bin
res.GCcounts = autoGCcounts("files.RData", "bins.RData", other.resources=list(account='rrg-bourqueg-ad'), genome=genome)

## QC (optional)
res.forQC = autoExtra("files.RData", "bins.RData", do=1, other.resources=list(account='rrg-bourqueg-ad')))
qc.samples.cluster(res.forQC) ## Run locally because it opens an interactive web browser application
##

## Normalize and call CNVs
res.df = autoNormTest("files.RData", "bins.RData", other.resources=list(account='rrg-bourqueg-ad')))
write.table(res.df, file='PopSV-CNVcalls.tsv', sep='\t', row.names=FALSE, quote=FALSE)

## Filter CNVs
res.filt.df = sv.summary.interactive(res.df) ## Run locally because it opens an interactive web browser application
write.table(res.filt.df, file='PopSV-CNVcalls-filtered.tsv', sep='\t', row.names=FALSE, quote=FALSE)



##
## Optional: Run additional samples using references from the previous analysis
##

## Option 1: in the same folder but using suffixes for the new batches
bam.files2 = read.table("bams2.tsv", as.is=TRUE, header=TRUE)
files.df = init.filenames(bam.files2, code="example2")
save(files.df, file="files2.RData")

res2.GCcounts = autoGCcorrect("files2.RData", "bins.RData", skip=1, file.suffix='batch2') # different suffix for batch2
res2.df = autoNormTest("files2.RData", "bins.RData", file.suffix.ref='', file.suffix='batch2') # and also specify suffix for reference analysis
write.table(res2.df, file='PopSV-CNVcalls-batch2.tsv', sep='\t', row.names=FALSE, quote=FALSE)

## Option 2: new batch in a separate folder
## Assuming that we work in a new "batch2" folder containing the "bams2.tsv' file
setwd('batch2') # update working directory
bam.files2 = read.table("bams2.tsv", as.is=TRUE, header=TRUE)
files.df = init.filenames(bam.files2, code="example2")
save(files.df, file="files2.RData")

res.GCcounts = autoGCcorrect("files2.RData", "../bins.RData", skip=1)
res.df = autoNormTest("files2.RData", "../bins.RData", ref.dir='..') # ref.dir specify the folder containing the reference analysis
write.table(res.df, file='PopSV-CNVcalls-batch2.tsv', sep='\t', row.names=FALSE, quote=FALSE)
