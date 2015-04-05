context("QC and merge reference samples")

## Simulate 3 samples
bin.df =  data.frame(chr=sample(1:22,1e2, replace=TRUE), start=seq(1,1e4,100), end=seq(99,9999,100))
samps = paste0("samp",1:3)
samps.f = paste(samps, "-bc.tsv")
sapply(1:3, function(ii){
  bin.df$bc = rnorm(100,1000,100)
  write.table(bin.df, file=samps.f[ii], row.names=FALSE, quote=FALSE, sep="\t")
})
files.df = data.frame(sample=samps, bc.gc.gz=comp.index.files(samps.f))

test_that("No missing columns",{})
test_that("Can pick a subset of reference samples",{})
test_that("Writes output files correctly",{})
test_that("Works with different chunk sizes",{})
test_that("Works with different columns",{})
test_that("Broad median normalization applied",{})
test_that("Samples and coordinates are correctly merged",{})
test_that("PCA is computed",{})
test_that("It's robust to NA/Infinite values",{})
