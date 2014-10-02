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
submitJobs(getGC.reg, 1, resources=list(walltime="30:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(getGC.reg)


