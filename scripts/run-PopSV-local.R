## This version is UNPRACTICAL AND NOT RECOMMENDED FOR REAL ANALYSIS. It is there as a more readable worflow. The BatchJobs version should be used in real analysis. Eventually, if BatchJobs is not desired, this script could help the user to deconstruct it to fit his computing cluster favorite strategy (but seriously don't waste your time and "learn" BatchJobs).

## FTE : FOR TOY EXAMPLE: lines marked with this flag should be adapted to the analysis. Here I used smaller numbers to run a small toy example. '=>' shows what should be used instead for an analysis in practice.

## Installation
## devtools::install_github("jmonlong/PopSV")

library(PopSV)

setwd("../../PopSV-toyex") ## FTE: set working directory path => put WD of the project
bam.files = read.table("bams.tsv", as.is=TRUE, header=TRUE)
bin.size = 1e3


## 1) Init file names and construct bins
files.df = init.filenames(bam.files, code="example")
bins.df = fragment.genome.hp19(bin.size)
bins.df = subset(bins.df, chr==22) ## FTE: chr 22 only => remove line
save(bins.df, file="bins.RData")
rm(bins.df)


## 2) Get GC content in each bin
bins.df = getGC.hg19(bins.df)


## 3-4) Get bin counts in each sample and GC correction
bc.o = sapply(1:nrow(files.df), function(file.i){
    bin.bam(files.df$bam[file.i], bins.df, files.df$bc[file.i])
    correct.GC(files.df$bc.gz[file.i], bins.df, files.df$bc.gc[file.i])
})


## 5) Sample QC and reference sample definition
## You may(should) have an idea of which samples can be used as reference (e.g. normal in a cancer project; controls in a case/control project; eventually every samples)
ref.samples = subset(files.df, group=="normal")$sample ## Here I had this information in the original bam info file in a 'group' column.

## If you have too many reference samples (lucky you) and want to find, say, 200 of them to use for the analysis, use 'nb.ref.samples=200' parameter in 'qc.samples(...)'

## OPTIONAL: If you suspect important batch effects that could create distinct groups of samples, the samples can be clustered first. Eventually analysis can be run separately on each cluster. THIS IS USUALLY NOT NECESSARY but is safe to check..
bc.rand = quick.count(subset(files.df, group=="normal"), bins.df, nb.cores=3, col.files="bc.gc.gz", nb.rand.bins=1e3) ## Gets counts for all samples on a subset of 1000 bins
qc.samples.cluster(bc.rand) ## Run locally because it opens an interactive web browser apllication
## END OPTIONAL

bc.all.f = "bc-gcCor-all.tsv"
pdf("sampQC.pdf")
samp.qc.o = qc.samples(files.df, bins.df, ref.samples, outfile.prefix=bc.all.f)
dev.off()


## 6) Normalize and compute Z scores

## 6a - Using targeted normalization (Recommended)

#### Normalize reference samples
bins.df = chunk.bin(bins.df, bg.chunk.size=2e4, sm.chunk.size=1e4, large.chr.chunks=TRUE) ## FTE: smaller chunks => 'bg.chunk.size=1e5' recommended
save(bins.df, file="bins.RData")

norm.o = lapply(unique(bins.df$sm.chunk), function(chunk.id){
  bc.df = read.bedix(samp.qc.o$bc, subset(bins.df, bg.chunk==subset(bins.df, sm.chunk==chunk.id)$bg.chunk[1]))
  tn.norm(bc.df, samp.qc.o$cont.sample, bins=subset(bins.df, sm.chunk==chunk.id)$bin)
})
out.files = paste("ref", c("bc-norm.tsv", "msd.tsv"), sep="-")
file.remove(out.files)
tmp = lapply(norm.o, function(res){
    write.table(res$bc.norm, file=out.files[1], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[1]), col.names=!file.exists(out.files[1]))
    write.table(res$norm.stats, file=out.files[2], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[2]), col.names=!file.exists(out.files[2]))
})

#### Compute Z-scores in reference samples
res.z = z.comp(out.files[1], ref.samples, out.files[2])
write.split.samples(res.z, files.df, ref.samples, res.n=c("z","fc"), files.col=c("z","fc"), compress.index=TRUE)


#### Normalization and Z-score computation for other samples
cases.o = lapply(setdiff(files.df$sample, ref.samples), function(samp){
  tn.test.sample(samp, files.df, samp.qc.o$cont.sample, samp.qc.o$bc, out.files[2])
})



## 7) Abnormal bin calling
call.o = lapply(files.df$sample, function(samp){
  call.abnormal.cov(files.df=files.df, samp=samp, out.pdf=paste0(samp,"/",samp,"-abCovCall.pdf"), FDR.th=.01, merge.cons.bins="stitch", z.th="sdest2N")
})

res.df = plyr::ldply(call.o, identity)

## Open locally because the output opens an interactive web browser application
sv.summary.interactive(res.df)
##
