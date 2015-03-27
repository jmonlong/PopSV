context("Tabix files manipulation")

bin.df =  data.frame(chr=sample(1:22,1e3, replace=TRUE), start=seq(1,1e5,100), end=seq(99,99999,100))
bin.df = bin.df[order(bin.df$chr, bin.df$start),]

b1.df = cbind(bin.df, bc=rnorm(nrow(bin.df)))
b2.df = cbind(bin.df, bc=rnorm(nrow(bin.df)), bc2=rnorm(nrow(bin.df)))

test_that("Writes and reads entire files",{
  write.table(b1.df, file="test.tsv", quote=FALSE, sep="\t", row.names=FALSE)
  write.table(b2.df, file="test2.tsv", quote=FALSE, sep="\t", row.names=FALSE)
  expect_equal(length(comp.index.files(c("test.tsv","test2.tsv"))), 2)
  df = read.bedix("test.tsv.bgz")
  expect_equal(df$bc, b1.df$bc)
  expect_equal(file.remove(c("test.tsv.bgz","test2.tsv.bgz","test.tsv.bgz.tbi","test2.tsv.bgz.tbi")), rep(TRUE,4))
})

test_that("Removes temporary files or not",{
  write.table(b1.df, file="test.tsv", quote=FALSE, sep="\t", row.names=FALSE)
  expect_equal(length(comp.index.files("test.tsv", rm.input=FALSE)), 1)
  expect_equal(file.exists("test.tsv"), TRUE)
  expect_equal(length(comp.index.files("test.tsv")), 1)
  expect_equal(file.exists("test.tsv"), FALSE)
  file.remove("test.tsv.bgz","test.tsv.bgz.tbi")
})

test_that("Reads subsets",{
  write.table(b1.df, file="test.tsv", quote=FALSE, sep="\t", row.names=FALSE)
  expect_equal(length(comp.index.files("test.tsv")), 1)
  ii = sort(sample.int(10))
  df = read.bedix("test.tsv.bgz", b1.df[ii,])
  expect_equal(df$bc, b1.df$bc[ii])
  ii = sort(sample.int(100))
  df = read.bedix("test.tsv.bgz", b1.df[ii,])
  expect_equal(df$bc, b1.df$bc[ii])
  ii = sort(sample.int(500))
  df = read.bedix("test.tsv.bgz", b1.df[ii,])
  expect_equal(df$bc, b1.df$bc[ii])
  file.remove("test.tsv.bgz","test.tsv.bgz.tbi")
})

test_that("Bins are ordered",{
  write.table(b1.df, file="test.tsv", quote=FALSE, sep="\t", row.names=FALSE)
  expect_equal(length(comp.index.files("test.tsv")), 1)
  ii = sample.int(100)
  df = read.bedix("test.tsv.bgz", b1.df[ii,])
  expect_equal(order(df$chr, df$start), 1:nrow(df))  
  file.remove("test.tsv.bgz","test.tsv.bgz.tbi")
})

test_that("Reads files with multiple columns",{
  write.table(b2.df, file="test.tsv", quote=FALSE, sep="\t", row.names=FALSE)
  expect_equal(length(comp.index.files("test.tsv")), 1)
  df = read.bedix("test.tsv.bgz")
  expect_equal(df$bc, b2.df$bc)
  expect_equal(df$bc2, b2.df$bc2)
  file.remove("test.tsv.bgz","test.tsv.bgz.tbi")
})

test_that("reads file with or without headers",{
  write.table(b2.df, file="test.tsv", quote=FALSE, sep="\t", row.names=FALSE)
  expect_equal(length(comp.index.files("test.tsv")), 1)
  expect_equal(nrow(read.bedix("test.tsv.bgz")), nrow(b2.df))
  write.table(b2.df, file="test.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  expect_equal(length(comp.index.files("test.tsv")), 1)
  expect_equal(nrow(read.bedix("test.tsv.bgz", header=FALSE)), nrow(b2.df))  
  file.remove("test.tsv.bgz","test.tsv.bgz.tbi")
})

test_that("reads file with or without .bgz extension",{})

test_that("Checks that files exists",{})

test_that("Checks that files are indexed",{})

test_that("Checks for bins missing columns",{})

