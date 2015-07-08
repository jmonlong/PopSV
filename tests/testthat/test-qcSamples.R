context("QC and merge reference samples")

## Simulate 3 samples
bin.df =  data.frame(chr=sample(1:22,1e2, replace=TRUE), start=seq(1,1e4,100), end=seq(99,9999,100))
samps = paste0("samp",1:4)
samps.f = paste0(samps, "-bc.tsv")
samp.med = c(1000,2000,3000,5000)
sapply(1:4, function(ii){
  bin.df = bin.df[order(bin.df$chr, bin.df$start),]
  bin.df$bc = round(rnorm(100,samp.med[ii],100),2)
  write.table(bin.df, file=samps.f[ii], row.names=FALSE, quote=FALSE, sep="\t")
})
files.df = data.frame(sample=samps, bc.gc.gz=comp.index.files(samps.f))

test_that("No missing columns",{
  expect_error(qc.samples(files.df, bin.df, col.bc="fakeCol"), "Missing column")
  files.bad = files.df
  files.bad$sample = NULL
  expect_error(qc.samples(files.bad, bin.df), "Missing column")
  bad.df = bin.df[,c("chr","start")]
  expect_error(qc.samples(files.df, bad.df), "Missing column")
  bad.df = bin.df[,c("end","start")]
  expect_error(qc.samples(files.df, bad.df), "Missing column")
  bad.df = bin.df[,c("chr","end")]
  expect_error(qc.samples(files.df, bad.df), "Missing column")
})

test_that("Can pick a subset of reference samples",{
  qc.o = qc.samples(files.df, bin.df, ref.samples = samps[1:3], outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE)
  expect_equal(qc.o$ref.samples, samps[1:3])
  expect_true(file.remove("temp.tsv"))
  qc.o = qc.samples(files.df, bin.df, nb.ref.samples = 3, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE)
  expect_equal(length(qc.o$ref.samples), 3)
  expect_true(file.remove("temp.tsv"))
})

test_that("Writes output files correctly",{
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE)
  bcm = read.table(qc.o$bc, header=TRUE, as.is=TRUE)
  expect_equal(dim(bcm), c(nrow(bin.df), 7))
  expect_true(file.remove("temp.tsv"))
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=TRUE)
  bin.sub = bin.df[sample.int(nrow(bin.df), 20),]
  bcm = read.bedix(qc.o$bc, bin.sub)
  expect_equal(dim(bcm), c(nrow(bin.sub), 7))
  expect_true(file.remove("temp.tsv.bgz"))
  expect_true(file.remove("temp.tsv.bgz.tbi"))
})


test_that("Works with different chunk sizes",{
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE, chunk.size = 30)
  bcm = read.table(qc.o$bc, header=TRUE, as.is=TRUE)
  expect_equal(dim(bcm), c(nrow(bin.df), 7))
  expect_true(file.remove("temp.tsv"))
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE, chunk.size = 50)
  bcm = read.table(qc.o$bc, header=TRUE, as.is=TRUE)
  expect_equal(dim(bcm), c(nrow(bin.df), 7))
  expect_true(file.remove("temp.tsv"))
})

test_that("Works with different columns",{
  files.df$otherCol = files.df$bc.gc.gz
  files.df$bc.gc.gz = NULL
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE, chunk.size = 50, col.bc="otherCol")
  bcm = read.table(qc.o$bc, header=TRUE, as.is=TRUE)
  expect_equal(dim(bcm), c(nrow(bin.df), 7))
  expect_true(file.remove("temp.tsv"))
})

test_that("Broad median normalization applied",{
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE, chunk.size = 20)
  bcm = read.table(qc.o$bc, header=TRUE, as.is=TRUE)
  bcm.med = apply(bcm[,-(1:3)], 2, median, na.rm=TRUE)
  expect_true(all(abs(bcm.med - median(bcm.med))<100))
  expect_true(file.remove("temp.tsv"))
})

test_that("Samples and coordinates are correctly merged",{
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE, chunk.size = 20)
  bcm = read.table(qc.o$bc, header=TRUE, as.is=TRUE)
  expect_true(all(samps %in% colnames(bcm)))
  bc2 = read.table(paste0(samps.f[2], ".bgz"), header=TRUE, as.is=TRUE)
  bc2 = bc2[order(bc2$chr, bc2$start), ]
  expect_equal(rank(bc2$bc), rank(bcm[,samps[2]]))
  expect_true(file.remove("temp.tsv"))
})

test_that("PCA is computed",{
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE, chunk.size = 20)
  expect_true(nrow(qc.o$pca.ref.df)>0)
  expect_null(qc.o$pca.all.df)
  expect_true(file.remove("temp.tsv"))
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", nb.ref.samples=3, plot=FALSE, appendIndex.outfile=FALSE, chunk.size = 20)
  expect_true(nrow(qc.o$pca.all.df)>0)
  expect_true(file.remove("temp.tsv"))
})

test_that("It's robust to NA/Infinite values",{
  ii = 3
  bin.df = bin.df[order(bin.df$chr, bin.df$start),]
  bin.df$bc = rnorm(100,samp.med[ii],100)
  bin.df$bc[1:2] = c(NA, Inf)
  write.table(bin.df, file=samps.f[ii], row.names=FALSE, quote=FALSE, sep="\t")
  comp.index.files(samps.f[ii])
  qc.o = qc.samples(files.df, bin.df, outfile.prefix="temp.tsv", plot=FALSE, appendIndex.outfile=FALSE, chunk.size = 20)
  expect_true(nrow(qc.o$pca.ref.df)>0)
  expect_true(file.remove("temp.tsv"))
})

## Clean
files = unique(files.df$bc.gc.gz)
file.remove(c(as.character(files), paste0(files, ".tbi")))
