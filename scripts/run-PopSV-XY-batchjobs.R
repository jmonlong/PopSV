## FTE : FOR TOY EXAMPLE: lines marked with this flag should be adapted to the analysis. Here I used smaller numbers to run a small toy example. '=>' shows what should be used instead for an real analysis.

## Installation
## devtools::install_github("jmonlong/PopSV", ref="forPaperXY")

library(BatchJobs)
library(PopSV)

setwd("../../PopSV-toyex") ## FTE: set working directory path => put working directory of the project
bam.files = read.table("bams.tsv", as.is=TRUE, header=TRUE)
bin.size = 1e3


## 1) Init file names and construct bins
files.df = init.filenames(bam.files, code="example")
bins.df = fragment.genome.hp19(bin.size, XY.chr=TRUE)
save(bins.df, file="bins.RData")
rm(bins.df)


## 2) Get GC content in each bin
## system("rm -rf getGC-files")
getGC.reg <- makeRegistry(id="getGC")
getGC.f <- function(imF){
    load(imF)
    library(PopSV)
    getGC.hg19(bins.df)
}
batchMap(getGC.reg, getGC.f,"bins.RData")
submitJobs(getGC.reg, findNotDone(getGC.reg), resources=list(walltime="1:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(getGC.reg)


## 3) Get bin counts in each sample and correct for GC bias
## system("rm -rf getBC-files")
getBC.reg <- makeRegistry(id="getBC")
getBC.f <- function(file.i, gc.reg, files.df){
  library(PopSV)
  bins.df = loadResult(gc.reg,1)
  bb.o = bin.bam(files.df$bam[file.i], bins.df, files.df$bc[file.i])
  correct.GC(files.df$bc.gz[file.i], bins.df, files.df$bc.gc[file.i])
  bb.o
}
batchMap(getBC.reg, getBC.f,1:nrow(files.df), more.args=list(gc.reg=getGC.reg, files.df=files.df))
submitJobs(getBC.reg, findNotDone(getBC.reg), resources=list(walltime="20:0:0", nodes="1", cores="1"), wait=function(retries) 100, max.retries=10)
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
bc.all.f = "bc-gcCor-all.tsv"
sampQC.pdf.f = "sampQC.pdf"
## system("rm -rf sampQC-files")
sampQC.reg <- makeRegistry(id="sampQC")
sampQC.f <- function(bc.all.f, bins.f, files.df, ref.samples, sampQC.pdf.f, nbc){
  library(PopSV)
  load(bins.f)
  pdf(sampQC.pdf.f)
  qc.o= qc.samples(files.df, bins.df, ref.samples, outfile.prefix=bc.all.f, nb.cores=nbc)
  dev.off()
  qc.o
}
batchMap(sampQC.reg, sampQC.f,bc.all.f, more.args=list(bins.f="bins.RData", files.df=files.df, ref.samples = ref.samples, sampQC.pdf.f=sampQC.pdf.f, nbc=3))
submitJobs(sampQC.reg, 1, resources=list(walltime="3:0:0", nodes="1", cores="3",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(sampQC.reg)
samp.qc.o = loadResult(sampQC.reg, 1)




######
##  In order to analyze chromosome X and Y, the analysis is split in 3 parts : all samples on autosomes, X on female samples, XY on male samples.
######
## The normalization step can be run simultaneously for the 3 parts. However, step 6 and 7 (z-scores and cases) should be run one at a time, starting by the automsome part (because the results of the X and XY parts will be appended to the files created in the first part).

### To be nice, this part could be run on an interactive node
bins.df = load("bins.RData")
bins.df = chunk.bin(bins.df, bg.chunk.size=2e4, sm.chunk.size=1e4, large.chr.chunks=TRUE) ## FTE: smaller chunks => 'bg.chunk.size=1e5' recommended
bins.sm.chunks = unique(bins.df$sm.chunk)
save(bins.df, file="bins.RData")
### End of the "nice guy interactive node" part

###
##  All samples on autosomes

## 5) Normalize and compute Z scores
## 5a - Using targeted normalization (Recommended)
#### Normalize reference samples
## system("rm -rf bcNormTN122-files")
bcNormTN122.reg <- makeRegistry(id="bcNormTN122")
bcNormTN122.f <- function(chunk.id, file.bc, file.bin, cont.sample){
    load(file.bin)
    library(PopSV)
    bins.df = subset(bins.df, chr %in% 1:22)
    bc.df = read.bedix(file.bc, subset(bins.df, bg.chunk==subset(bins.df, sm.chunk==chunk.id)$bg.chunk[1]))
    tn.norm(bc.df, cont.sample, bins=subset(bins.df, sm.chunk==chunk.id)$bin)
}
batchMap(bcNormTN122.reg, bcNormTN122.f,bins.sm.chunks, more.args=list(file.bc=samp.qc.o$bc, file.bin="bins.RData",cont.sample=samp.qc.o$cont.sample))
submitJobs(bcNormTN122.reg, findJobs(bcNormTN122.reg) , resources=list(walltime="10:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(bcNormTN122.reg)

## Write normalized bin counts and reference metrics
out.files.122 = paste("ref-1-22", c("bc-norm.tsv", "msd.tsv"), sep="-")
file.remove(out.files.122)
tmp = reduceResultsList(bcNormTN122.reg, fun=function(res, job){
    write.table(res$bc.norm, file=out.files.122[1], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files.122[1]), col.names=!file.exists(out.files.122[1]))
    write.table(res$norm.stats, file=out.files.122[2], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files.122[2]), col.names=!file.exists(out.files.122[2]))
})

## 6) Compute Z-scores in reference samples
## system("rm -rf zRef122-files")
zRef122.reg <- makeRegistry(id="zRef122")
zRef122.f <- function(bc.f, files.df){
  library(PopSV)
  z.comp(bc.f=bc.f, files.df=files.df, nb.cores=3, z.poisson=TRUE, chunk.size=1e4)
}
batchMap(zRef122.reg, zRef122.f,out.files.122[1], more.args=list(files.df=files.df))
submitJobs(zRef122.reg, 1, resources=list(walltime="6:0:0", nodes="1", cores="3",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(zRef122.reg)

## 7) Normalization and Z-score computation for other samples
## system("rm -rf callsCases122-files")
callsCases122.reg <- makeRegistry(id="callsCases122")
callsCases122.f <- function(samp, cont.sample, files.df, norm.stats.f, bc.ref.f){
  library(PopSV)
  tn.test.sample(samp, files.df, cont.sample, bc.ref.f, norm.stats.f, z.poisson=TRUE, compress.index=FALSE)
}
batchMap(callsCases122.reg, callsCases122.f,setdiff(files.df$sample, ref.samples), more.args=list(cont.sample=samp.qc.o$cont.sample, files.df=files.df, norm.stats.f=out.files.122[2], bc.ref.f=samp.qc.o$bc))
submitJobs(callsCases122.reg, findNotDone(callsCases122.reg), resources=list(walltime="20:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(callsCases122.reg)


###
##  chr X on females
ref.samples.female = subset(files.df, group=="normal" & gender=="female")$sample
cont.sample.female = sample(ref.samples.female, 1)

## 5) Normalize and compute Z scores
## 5a - Using targeted normalization (Recommended)
#### Normalize reference samples
## system("rm -rf bcNormTNXF-files")
bcNormTNXF.reg <- makeRegistry(id="bcNormTNXF")
bcNormTNXF.f <- function(chunk.id, file.bc, file.bin, samples, cont.sample){
    load(file.bin)
    library(PopSV)
    if(!any(bins.df$sm.chunk==chunk.id & bins.df$chr=="X")) return(NULL)
    bc.df = read.bedix(file.bc, subset(bins.df, bg.chunk==subset(bins.df, sm.chunk==chunk.id)$bg.chunk[1]))
    bc.df = bc.df[, c("chr","start","end", intersect(samples, colnames(bc.df)))]
    tn.norm(bc.df, cont.sample, bins=subset(bins.df, sm.chunk==chunk.id & chr=="X")$bin)
}
batchMap(bcNormTNXF.reg, bcNormTNXF.f,bins.sm.chunks, more.args=list(file.bc=samp.qc.o$bc, file.bin="bins.RData",samples=ref.samples.female, cont.sample=cont.sample.female))
submitJobs(bcNormTNXF.reg, findJobs(bcNormTNXF.reg) , resources=list(walltime="10:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(bcNormTNXF.reg)

## Write normalized bin counts and reference metrics
out.files.XF = paste("ref-X-female", c("bc-norm.tsv", "msd.tsv"), sep="-")
file.remove(out.files.XF)
tmp = reduceResultsList(bcNormTNXF.reg, fun=function(res, job){
    write.table(res$bc.norm, file=out.files.XF[1], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files.XF[1]), col.names=!file.exists(out.files.XF[1]))
    write.table(res$norm.stats, file=out.files.XF[2], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files.XF[2]), col.names=!file.exists(out.files.XF[2]))
})

## 6) Compute Z-scores in reference samples
## system("rm -rf zRefXF-files")
zRefXF.reg <- makeRegistry(id="zRefXF")
zRefXF.f <- function(bc.f, files.df){
  library(PopSV)
  z.comp(bc.f=bc.f, files.df=files.df, nb.cores=3, z.poisson=TRUE, chunk.size=1e4, out.msd.f=NULL, append=TRUE)
}
batchMap(zRefXF.reg, zRefXF.f,out.files.XF[1], more.args=list(files.df=files.df))
submitJobs(zRefXF.reg, 1, resources=list(walltime="6:0:0", nodes="1", cores="3",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(zRefXF.reg)

## 7) Normalization and Z-score computation for other samples
## system("rm -rf callsCasesXF-files")
callsCasesXF.reg <- makeRegistry(id="callsCasesXF")
callsCasesXF.f <- function(samp, cont.sample, files.df, norm.stats.f, bc.ref.f){
  library(PopSV)
  tn.test.sample(samp, files.df, cont.sample, bc.ref.f, norm.stats.f, z.poisson=TRUE, compress.index=FALSE, append=TRUE)
}
batchMap(callsCasesXF.reg, callsCasesXF.f,subset(files.df, gender=="female" & group!="normal")$sample, more.args=list(cont.sample=cont.sample.female, files.df=files.df, norm.stats.f=out.files.XF[2], bc.ref.f=samp.qc.o$bc))
submitJobs(callsCasesXF.reg, findNotDone(callsCasesXF.reg), resources=list(walltime="20:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(callsCasesXF.reg)


###
##  chr XY on males
ref.samples.male = subset(files.df, group=="normal" & gender=="male")$sample
cont.sample.male = sample(ref.samples.male, 1)

## 5) Normalize and compute Z scores
## 5a - Using targeted normalization (Recommended)
#### Normalize reference samples
## system("rm -rf bcNormTNXYM-files")
bcNormTNXYM.reg <- makeRegistry(id="bcNormTNXYM")
bcNormTNXYM.f <- function(chunk.id, file.bc, file.bin, samples, cont.sample){
    load(file.bin)
    library(PopSV)
    if(!any(bins.df$sm.chunk==chunk.id & bins.df$chr %in% c("X","Y"))) return(NULL)
    bc.df = read.bedix(file.bc, subset(bins.df, bg.chunk==subset(bins.df, sm.chunk==chunk.id)$bg.chunk[1]))
    bc.df = bc.df[, c("chr","start","end", intersect(samples, colnames(bc.df)))]
    tn.norm(bc.df, cont.sample, bins=subset(bins.df, sm.chunk==chunk.id & chr %in% c("X","Y"))$bin)
}
batchMap(bcNormTNXYM.reg, bcNormTNXYM.f,bins.sm.chunks, more.args=list(file.bc=samp.qc.o$bc, file.bin="bins.RData",samples=ref.samples.male, cont.sample=cont.sample.male))
submitJobs(bcNormTNXYM.reg, findJobs(bcNormTNXYM.reg) , resources=list(walltime="10:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(bcNormTNXYM.reg)

## Write normalized bin counts and reference metrics
out.files.XYM = paste("ref-X-male", c("bc-norm.tsv", "msd.tsv"), sep="-")
file.remove(out.files.XYM)
tmp = reduceResultsList(bcNormTNXYM.reg, fun=function(res, job){
    write.table(res$bc.norm, file=out.files.XYM[1], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files.XYM[1]), col.names=!file.exists(out.files.XYM[1]))
    write.table(res$norm.stats, file=out.files.XYM[2], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files.XYM[2]), col.names=!file.exists(out.files.XYM[2]))
})

## 6) Compute Z-scores in reference samples
## system("rm -rf zRefXYM-files")
zRefXYM.reg <- makeRegistry(id="zRefXYM")
zRefXYM.f <- function(bc.f, files.df){
  library(PopSV)
  z.comp(bc.f=bc.f, files.df=files.df, nb.cores=3, z.poisson=TRUE, chunk.size=1e4, out.msd.f=NULL, append=TRUE, compress.index=FALSE)
}
batchMap(zRefXYM.reg, zRefXYM.f,out.files.XYM[1], more.args=list(files.df=files.df))
submitJobs(zRefXYM.reg, 1, resources=list(walltime="6:0:0", nodes="1", cores="3",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(zRefXYM.reg)

## 7) Normalization and Z-score computation for other samples
## system("rm -rf callsCasesXYM-files")
callsCasesXYM.reg <- makeRegistry(id="callsCasesXYM")
callsCasesXYM.f <- function(samp, cont.sample, files.df, norm.stats.f, bc.ref.f){
  library(PopSV)
  tn.test.sample(samp, files.df, cont.sample, bc.ref.f, norm.stats.f, z.poisson=TRUE, append=TRUE)
}
batchMap(callsCasesXYM.reg, callsCasesXYM.f,subset(files.df, gender=="male" & group!="normal")$sample, more.args=list(cont.sample=cont.sample.male, files.df=files.df, norm.stats.f=out.files.XYM[2], bc.ref.f=samp.qc.o$bc))
submitJobs(callsCasesXYM.reg, findNotDone(callsCasesXYM.reg), resources=list(walltime="20:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(callsCasesXYM.reg)



#######
##   For all samples abnormal regions are called
#######

## 8) Abnormal bin calling
## system("rm -rf abCovCallCases-files")
abCovCallCases.reg <- makeRegistry(id="abCovCallCases")
abCovCallCases.f <- function(samp, files.df, norm.stats.f, bin.size){
  library(PopSV)
  call.abnormal.cov(files.df=files.df, samp=samp, out.pdf=paste0(samp,"/",samp,"-sdest2N-abCovCall.pdf"), FDR.th=.01, merge.cons.bins="stitch", z.th="sdest2N", norm.stats=norm.stats.f, stitch.dist=bin.size*2+1)
}
batchMap(abCovCallCases.reg, abCovCallCases.f, files.df$sample, more.args=list(files.df=files.df, norm.stats.f=out.files.122[2], bin.size=bin.size))
submitJobs(abCovCallCases.reg, findNotDone(abCovCallCases.reg) , resources=list(walltime="2:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(abCovCallCases.reg)

res.df = plyr::ldply(reduceResultsList(abCovCallCases.reg), identity)
res.df$.id = NULL
save(res.df, file="cnvs-FDR01-mergeStitch-thSdest2N.RData")

## Open locally because the output opens an interactive web browser application
sv.summary.interactive(res.df)
##
