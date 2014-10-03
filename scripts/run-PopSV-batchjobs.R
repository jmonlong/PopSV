library(BatchJobs)
library(PopSV)

bam.files = NULL
bins.df = NULL

## 1) Init file names
files.df = initFileNames(bam.files, code="example")

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
    bin.df = loadResult(gc.reg,1)
    binBam(files.df$bam[file.i], bin.df, files.df$bc[file.i])
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
    binBam(files.df$bc.gz[file.i], gc.df, files.df$bc.gc[file.i])
}
batchMap(gcCor.reg, gcCor.f,1:nrow(files.df), more.args=list(gc.reg=getGC.reg, files.df=files.df))
submitJobs(gcCor.reg, findNotDone(gcCor.reg), resources=list(walltime="12:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(gcCor.reg)
all(reduceResultsVector(gcCor.reg)==files.df$bc.gc.gz)
