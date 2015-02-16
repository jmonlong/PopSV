## FTE : FOR TOY EXAMPLE: lines marked with this flag should be adapted to the analysis. Here I used smaller numbers to run a small toy example. '=>' shows what should be used instead for an analysis in practice.

library(BatchJobs)
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
## system("rm -rf getGC-files")
getGC.reg <- makeRegistry(id="getGC")
getGC.f <- function(imF){
    load(imF)
    library(PopSV)
    getGC.hg19(bins.df)
}
batchMap(getGC.reg, getGC.f,"bins.RData")
submitJobs(getGC.reg, 1, resources=list(walltime="3:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(getGC.reg)


## 3) Get bin counts in each sample
## system("rm -rf getBC-files")
getBC.reg <- makeRegistry(id="getBC")
getBC.f <- function(file.i, gc.reg, files.df){
    library(PopSV)
    bins.df = loadResult(gc.reg,1)
    bin.bam(files.df$bam[file.i], bins.df, files.df$bc[file.i])
}
batchMap(getBC.reg, getBC.f,1:nrow(files.df), more.args=list(gc.reg=getGC.reg, files.df=files.df))
submitJobs(getBC.reg, findNotDone(getBC.reg), resources=list(walltime="10:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(getBC.reg)

## OPTIONAL QC: check the total number of reads counted in all samples
pdf("qc-nbReads.pdf")
qplot(reduceResultsVector(getBC.reg, fun=function(job, res)res$nb.reads)) + xlab("total number of reads") + ylab("number of samples") + theme_bw()
dev.off()
##


## 4) GC correction
## TODO: Merge this step with the previous one !
## system("rm -rf gcCor")
gcCor.reg <- makeRegistry(id="gcCor")
gcCor.f <- function(file.i, gc.reg, files.df){
    library(PopSV)
    gc.df = loadResult(gc.reg,1)
    correct.GC(files.df$bc.gz[file.i], gc.df, files.df$bc.gc[file.i])
}
batchMap(gcCor.reg, gcCor.f,1:nrow(files.df), more.args=list(gc.reg=getGC.reg, files.df=files.df))
submitJobs(gcCor.reg, findNotDone(gcCor.reg), resources=list(walltime="3:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(gcCor.reg)


## 5) Sample QC and reference sample definition
## You may(should) have an idea of which samples can be used as reference (e.g. normal in a cancer project; controls in a case/control project; eventually every samples)
ref.samples = subset(files.df, group=="normal")$sample ## Here I had this information in the original bam info file in a 'group' column.
## If you have too many reference samples (lucky you) and want to find, say, 200 of them to use for the analysis, use 'nb.ref.samples=200' parameter in 'qc.samples(...)'
bc.all.f = "bc-gcCor-all.tsv"
sampQC.pdf.f = "sampQC.pdf"
## system("rm -rf sampQC-files")
sampQC.reg <- makeRegistry(id="sampQC")
sampQC.f <- function(bc.all.f, gc.reg, files.df, ref.samples, sampQC.pdf.f){
    library(PopSV)
    gc.df = loadResult(gc.reg,1)
    pdf(sampQC.pdf.f)
    qc.o= qc.samples(files.df, gc.df, ref.samples, outfile.prefix=bc.all.f, nb.cores=6)
    dev.off()
    qc.o
}
batchMap(sampQC.reg, sampQC.f,bc.all.f, more.args=list(gc.reg=getGC.reg, files.df=files.df, ref.samples = ref.samples, sampQC.pdf.f=sampQC.pdf.f))
submitJobs(sampQC.reg, 1, resources=list(walltime="3:0:0", nodes="1", cores="6",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(sampQC.reg)
samp.qc.o = loadResult(sampQC.reg, 1)


## 6) Normalize and compute Z scores

## 6a - Using targeted normalization (Recommended)

#### Normalize reference samples
## system("rm -rf bcNormTN-files")
bcNormTN.reg <- makeRegistry(id="bcNormTN")
### To be nice, this part could be run on an interactive node
bins.df = loadResult(getGC.reg,1)
bins.df = chunk.bin(bins.df, bg.chunk.size=2e4, sm.chunk.size=1e3, large.chr.chunks=TRUE) ## FTE: smaller chunks => 'bg.chunk.size=1e5' recommended
save(bins.df, file="bins.RData")
bcNormTN.f <- function(chunk.id, file.bc, file.bin, cont.sample){
    load(file.bin)
    library(PopSV)
    bc.df = read.bedix(file.bc, subset(bins.df, bg.chunk==subset(bins.df, sm.chunk==chunk.id)$bg.chunk[1]))
    tn.norm(bc.df, cont.sample, bins=subset(bins.df, sm.chunk==chunk.id)$bin)
}
batchMap(bcNormTN.reg, bcNormTN.f,unique(bins.df$sm.chunk), more.args=list(file.bc=samp.qc.o$bc, file.bin="bins.RData",cont.sample=samp.qc.o$cont.sample))
### End of the "nice guy interactive node" part
submitJobs(bcNormTN.reg, findJobs(bcNormTN.reg) , resources=list(walltime="12:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(bcNormTN.reg)
out.files = paste("ref", c("bc-norm.tsv", "msd.tsv"), sep="-")
file.remove(out.files)
tmp = reduceResultsList(bcNormTN.reg, fun=function(res, job){
    write.table(res$bc.norm, file=out.files[1], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[1]), col.names=!file.exists(out.files[1]))
    write.table(res$norm.stats, file=out.files[2], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[2]), col.names=!file.exists(out.files[2]))
})

#### Compute Z-scores in reference samples
## system("rm -rf zRef-files")
zRef.reg <- makeRegistry(id="zRef")
zRef.f <- function(file.bc, samples, msd.f, files.df){
    library(PopSV)
    res = z.comp(file.bc, samples, msd.f, nb.cores=6)
    write.split.samples(res, files.df, samples, res.n=c("z","fc"), files.col=c("z","fc"), compress.index=TRUE)
}
batchMap(zRef.reg, zRef.f,out.files[1], more.args=list(samples=ref.samples, msd.f=out.files[2], files.df=files.df))
submitJobs(zRef.reg, 1, resources=list(walltime="12:0:0", nodes="1", cores="6",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(zRef.reg)


#### Normalization and Z-score computation for other samples
