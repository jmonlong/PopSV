## FTE : FOR TOY EXAMPLE: lines marked with this flag should be adapted to the analysis. Here I used smaller numbers to run a small toy example. '=>' shows what should be used instead for an real analysis.

#### Installation
## devtools::install_github("jmonlong/PopSV")

library(BatchJobs)
library(PopSV)

setwd("../../PopSV-toyex") ## FTE: set working directory path => put working directory of the project
bam.files = read.table("bams.tsv", as.is=TRUE, header=TRUE)
bin.size = 1e3


## 1) Init file names and construct bins
files.df = init.filenames(bam.files, code="example")
bins.df = fragment.genome.hp19(bin.size)
bins.df = subset(bins.df, chr==22) ## FTE: chr 22 only => remove line
## It's good to save these files, they'll be used a lot, you can 'load()' them if you need them later
save(files.df, file="files.RData")
save(bins.df, file="bins.RData")


## 2) Get GC content in each bin
getGC.reg <- makeRegistry(id="getGC")
getGC.f <- function(imF){
  load(imF)
  library(PopSV)
  bins.df = getGC.hg19(bins.df)
  save(bins.df, file=imF)
}
batchMap(getGC.reg, getGC.f,"bins.RData")
submitJobs(getGC.reg, 1, resources=list(walltime="2:0:0", nodes="1", cores="1"))
showStatus(getGC.reg)


## 3) Get bin counts in each sample and correct for GC bias
getBC.reg <- makeRegistry(id="getBC")
getBC.f <- function(file.i, bins.f, files.df){
  library(PopSV)
  load(bins.f)
  bb.o = bin.bam(files.df$bam[file.i], bins.df, files.df$bc[file.i])
  correct.GC(files.df$bc.gz[file.i], bins.df, files.df$bc.gc[file.i])
  bb.o
}
batchMap(getBC.reg, getBC.f,1:nrow(files.df), more.args=list(bins.f="bins.RData", files.df=files.df))
submitJobs(getBC.reg, findNotDone(getBC.reg), resources=list(walltime="20:0:0", nodes="1", cores="1"))
showStatus(getBC.reg)

## OPTIONAL QC: check the total number of reads counted in all samples
library(ggplot2)
pdf("qc-nbReads.pdf")
qplot(reduceResultsVector(getBC.reg, fun=function(job, res)res$nb.reads)) + xlab("total number of reads") + ylab("number of samples") + theme_bw()
dev.off()
##
## OPTIONAL QC : If you suspect important batch effects that could create distinct groups of samples, the samples can be clustered first. Eventually analysis can be run separately on each cluster. THIS IS USUALLY NOT NECESSARY but is safe to check..
bc.rand = quick.count(subset(files.df, group=="normal"), bins.df, nb.cores=3, col.files="bc.gc.gz", nb.rand.bins=1e3) ## Gets counts for all samples on a subset of 1000 bins
qc.samples.cluster(bc.rand) ## Run locally because it opens an interactive web browser apllication
## END OPTIONAL

## 4) Sample QC and reference sample definition
## You may(should) have an idea of which samples can be used as reference (e.g. normal in a cancer project; controls in a case/control project; eventually every samples)
ref.samples = subset(files.df, group=="normal")$sample ## Here I had this information in the original bam info file in a 'group' column.
## If you have too many reference samples (lucky you) and want to find, say, 200 of them to use for the analysis, use 'nb.ref.samples=200' parameter in 'qc.samples(...)'
sampQC.reg <- makeRegistry(id="sampQC")
sampQC.f <- function(bins.f, files.df, ref.samples){
  library(PopSV)
  load(bins.f)
  pdf("sampQC.pdf")
  qc.o= qc.samples(files.df, bins.df, ref.samples, outfile.prefix="bc-gcCor-all.tsv", nb.cores=3)
  dev.off()
  qc.o
}
batchMap(sampQC.reg, sampQC.f,"bins.RData", more.args=list(files.df=files.df, ref.samples = ref.samples))
submitJobs(sampQC.reg, 1, resources=list(walltime="3:0:0", nodes="1", cores="3"))
showStatus(sampQC.reg)
samp.qc.o = loadResult(sampQC.reg, 1)


## 5) Normalize bin counts in reference samples
bcNormTN.reg <- makeRegistry(id="bcNormTN")
load("bins.RData")
bins.df = chunk.bin(bins.df, bg.chunk.size=2e4, sm.chunk.size=1e4, large.chr.chunks=TRUE) ## FTE: smaller chunks => 'bg.chunk.size=1e5' recommended
save(bins.df, file="bins.RData")
bcNormTN.f <- function(chunk.id, file.bc, file.bin, cont.sample){
  load(file.bin)
  library(PopSV)
  bc.df = read.bedix(file.bc, subset(bins.df, bg.chunk==subset(bins.df, sm.chunk==chunk.id)$bg.chunk[1]))
  tn.norm(bc.df, cont.sample, bins=subset(bins.df, sm.chunk==chunk.id)$bin)
}
batchMap(bcNormTN.reg, bcNormTN.f,unique(bins.df$sm.chunk), more.args=list(file.bc=samp.qc.o$bc, file.bin="bins.RData",cont.sample=samp.qc.o$cont.sample))
submitJobs(bcNormTN.reg, findExpired(bcNormTN.reg) , resources=list(walltime="12:0:0", nodes="1", cores="1"))
showStatus(bcNormTN.reg)

#### Write normalized bin counts and reference metrics
out.files = paste("ref", c("bc-norm.tsv", "msd.tsv"), sep="-")
file.remove(out.files)
tmp = reduceResultsList(bcNormTN.reg, fun=function(res, job){
    write.table(res$bc.norm, file=out.files[1], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[1]), col.names=!file.exists(out.files[1]))
    write.table(res$norm.stats, file=out.files[2], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[2]), col.names=!file.exists(out.files[2]))
})


## 6) Compute Z-scores in reference samples
zRef.reg <- makeRegistry(id="zRef")
zRef.f <- function(bc.f, files.df){
  library(PopSV)
  z.comp(bc.f=bc.f, files.df=files.df, nb.cores=3, z.poisson=TRUE, chunk.size=1e4)
}
batchMap(zRef.reg, zRef.f,out.files[1], more.args=list(files.df=files.df))
submitJobs(zRef.reg, 1, resources=list(walltime="6:0:0", nodes="1", cores="3"))
showStatus(zRef.reg)


## 7) Normalization and Z-score computation for other samples
normZcases.reg <- makeRegistry(id="normZcases")
normZcases.f <- function(samp, cont.sample, files.df, norm.stats.f, bc.ref.f){
  library(PopSV)
  tn.test.sample(samp, files.df, cont.sample, bc.ref.f, norm.stats.f, z.poisson=TRUE, aberrant.cases=TRUE)
}
batchMap(normZcases.reg, normZcases.f,setdiff(files.df$sample, ref.samples), more.args=list(cont.sample=samp.qc.o$cont.sample, files.df=files.df, norm.stats.f=out.files[2], bc.ref.f=samp.qc.o$bc))
submitJobs(normZcases.reg, findNotDone(normZcases.reg), resources=list(walltime="20:0:0", nodes="1", cores="1"))
showStatus(normZcases.reg)


## 8) Abnormal bin calling
abCovCallCases.reg <- makeRegistry(id="abCovCallCases")
abCovCallCases.f <- function(samp, files.df, norm.stats.f, bins.f, bin.size){
  library(PopSV)
  load(bins.f)
  call.abnormal.cov(files.df=files.df, samp=samp, out.pdf=paste0(samp,"-sdest-abCovCall.pdf"), FDR.th=.001, merge.cons.bins="stitch", z.th="sdest", norm.stats=norm.stats.f, stitch.dist=bin.size*2+1, gc.df = bins.df)
}
batchMap(abCovCallCases.reg, abCovCallCases.f, files.df$sample, more.args=list(files.df=files.df, norm.stats.f=out.files[2], bins.f="bins.RData",bin.size=bin.size))
submitJobs(abCovCallCases.reg, findNotDone(abCovCallCases.reg) , resources=list(walltime="1:0:0", nodes="1", cores="1"))
showStatus(abCovCallCases.reg)

res.df = do.call(rbind, reduceResultsList(abCovCallCases.reg), identity)
save(res.df, file="cnvs-FDR001-mergeStitch-thSdest.RData")

## Open locally because the output opens an interactive web browser application
sv.summary.interactive(res.df)
##
