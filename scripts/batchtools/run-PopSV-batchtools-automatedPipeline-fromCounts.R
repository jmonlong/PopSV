## If the bin counts are already available, you should make sure that the files
## are bgzipped and indexed. Then the bin counts should be GC corrected.

library(PopSV)
source('automatedPipeline-batchtools.R')

genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
## or
## genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

##
## Preparation
## Run only once to create the files files.RData and bins.RData
##

## Assuming a file with a column 'sample' with sample names and
## 'bincount' with path to the bin count file (columns chr/start/end/bc)
bc.files = read.table("bincounts.tsv", as.is=TRUE, header=TRUE)
bin.size = 1e3

#### Init file names and construct bins
bc.files$bam = NA
files.df = init.filenames(bc.files, code="example")
save(files.df, file="files.RData")
## If the bins are already counted, you should have the bins definition, e.g.
bins.df = read.table('bins.bed', sep='\t', as.is=TRUE)
colnames(bins.df) = c('chr', 'start', 'end')
save(bins.df, file="bins.RData")
####

## Order, bgzip and index
tmp = comp.index.files(files.df$bincount, files.df$bc, rm.input=FALSE, reorder=TRUE)
## Note: if reorder=TRUE it might be worth doing this in a job rather than
## the login node. If not, it's just doing bgzip and indexing.


##
## Analysis
## Can be stopped and restarted. No need to rerun the preparation commands
##

## Bin and count reads in each bin
res.GCcounts = autoGCcorrect("files.RData", "bins.RData", other.resources=list(account='rrg-bourqueg-ad'), genome=genome)

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

