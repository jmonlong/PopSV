library(batchtools)

## To start again from scratch, manually remove the 'test' folder.
reg <- makeRegistry('test', seed=123)
## reg = loadRegistry('test', writeable=TRUE) ## If the registry has already been created before

test.f <- function(ii){
  return(mean(rnorm(10,ii)))
}

batchMap(reg=reg, test.f, 1:2)
submitJobs(reg=reg, ids=findJobs(reg=reg), resources=list(walltime='10:00', cores=1))
waitForJobs(reg=reg, sleep=10)
getStatus(reg=reg)
reduceResultsList(reg=reg)
