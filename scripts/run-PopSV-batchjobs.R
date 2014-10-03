library(BatchJobs)
library(PopSV)
library(dplyr)

bam.files = NULL ## read.table("bamfiles.tsv", as.is=TRUE)
bin.size = 1e3

## 1) Init file names
files.df = initFileNames(bam.files, code="example")
bins.df = fragment.genome.hp19(bin.size)

## 2) Get GC content
## system("rm -rf getGC")
getGC.reg <- makeRegistry(id="getGC", seed=123, file.dir="getGC")
batchMap(getGC.reg, getGC.hg19,bins.df)
submitJobs(getGC.reg, 1, resources=list(walltime="12:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(getGC.reg)

## 3) Get bin counts
## system("rm -rf getBC")
getBC.reg <- makeRegistry(id="getBC", seed=123, file.dir="getBC")
getBC.f <- function(file.i, gc.reg, files.df){
    bins.df = loadResult(gc.reg,1)
    binBam(files.df$bam[file.i], bins.df, files.df$bc[file.i])
}
batchMap(getBC.reg, getBC.f,1:nrow(files.df), more.args=list(gc.reg=getGC.reg, files.df=files.df))
submitJobs(getBC.reg, findNotDone(getBC.reg), resources=list(walltime="30:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(getBC.reg)
all(reduceResultsVector(getBC.reg)==files.df$bc.gz)

## 4) GC correction
## system("rm -rf gcCor")
gcCor.reg <- makeRegistry(id="gcCor", seed=123, file.dir="gcCor")
gcCor.f <- function(file.i, gc.reg, files.df){
    gc.df = loadResult(gc.reg,1)
    correct.GC(files.df$bc.gz[file.i], gc.df, files.df$bc.gc[file.i])
}
batchMap(gcCor.reg, gcCor.f,1:nrow(files.df), more.args=list(gc.reg=getGC.reg, files.df=files.df))
submitJobs(gcCor.reg, findNotDone(gcCor.reg), resources=list(walltime="12:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(gcCor.reg)
all(reduceResultsVector(gcCor.reg)==files.df$bc.gc.gz)

## 5) Sample QC
### First guess at we the reference will be.
### (e.g. normal in  cancer study, controls in case/control study)
ref.samples = NULL
bc.all.f = "bc-gcCor-all.tsv"
sampQC.pdf.f = "sampQC.pdf"
## system("rm -rf sampQC")
sampQC.reg <- makeRegistry(id="sampQC", seed=123, file.dir="sampQC")
sampQC.f <- function(bc.all.f, gc.reg, files.df, ref.samples, sampQC.pdf.f){
    gc.df = loadResult(gc.reg,1)
    qc.samples(files.df, gc.df, ref.samples, bc.all.f, sampQC.pdf.f)
}
batchMap(sampQC.reg, sampQC.f,bc.all.f, more.args=list(gc.reg=getGC.reg,
                                            files.df=files.df, ref.samples = ref.samples,
                                            sampQC.pdf.f=sampQC.pdf.f))
submitJobs(sampQC.reg, 1, resources=list(walltime="12:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(sampQC.reg)
samp.qc.o = loadResult(sampQC.reg, 1)

## Potentially refine the reference samples definition
ref.samples = subset(samp.qc.o$dstat, Dstat > .8)$sample
cont.sample = arrange(samp.qc.o$dstat, desc(Dstat))$sample[1]

## 6) Normalize and compute Z scores
## system("rm -rf bcNormZ")
bcNormZ.reg <- makeRegistry(id="bcNormZ", seed=123, file.dir="bcNormZ")
chunk.size = 1e5
### To be nice, this part could be run on an interactive node
bins.df = loadResult(getGC.reg,1)
bins.df$chunk = sample(rep(1:ceiling(nrow(bins.df)/chunk.size),each=chunk.size)[1:nrow(bins.df)])
bcNormZ.f <- function(chunk.id, file.bc, imF){
    load(imF)
    bc.df = read.bedix(file.bc, subset(bins.df, chunk==chunk.id))
    tn.norm(bc.df, cont.sample, ref.samples)
}
imF = "bcNorm-temp.RData"
save(bins.df, cont.sample, ref.samples, file=imF)
batchMap(bcNormZ.reg, bcNormZ.f,1:max(bins.df$chunk), more.args=list(file.bc=samp.qc.o$bc, imF=imF))
### End of the "nice guy interactive node" part
submitJobs(bcNormZ.reg, findJobs(bcNormZ.reg) , resources=list(walltime="12:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(bcNormZ.reg)
