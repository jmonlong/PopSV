context("Bin counting from a BAM file")

bin.df = data.frame(chr="chr18", start=seq(1,15e5,1e4), stringsAsFactors=FALSE)
bin.df$end = bin.df$start + 1e4 -1
bam.f = "../dataForTests/NA12878_S1_chr18_1_1000000_s01.bam"

test_that("Arguments are consistent", {
  expect_error(bin.bam("fakefile.bam", bin.df, outfile.prefix=NULL, appendIndex.outfile=TRUE), "please provide 'outfile.prefix'")
  bad.df = bin.df[,c("chr","start")]
  expect_error(bin.bam(bam.f, bad.df, outfile.prefix=NULL, appendIndex.outfile=FALSE), "Missing column")
  bad.df = bin.df[,c("end","start")]
  expect_error(bin.bam(bam.f, bad.df, outfile.prefix=NULL, appendIndex.outfile=FALSE), "Missing column")
  bad.df = bin.df[,c("chr","end")]
  expect_error(bin.bam(bam.f, bad.df, outfile.prefix=NULL, appendIndex.outfile=FALSE), "Missing column")
  good.df = bin.df[,c("chr","end", "start")]
  expect_equal(length(bin.bam(bam.f, good.df, outfile.prefix=NULL, appendIndex.outfile=FALSE)), 2)
})

test_that("Files exists",{
  expect_error(bin.bam("fakefile.bam", bin.df, appendIndex.outfile=FALSE), "file not found")
  write("test", file="fakefile.bam")
  expect_error(bin.bam("fakefile.bam", bin.df, appendIndex.outfile=FALSE), "Index file is missing")
  file.remove("fakefile.bam")
})

test_that("Works with both chr naming",{
  expect_equal(length(bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE)), 2)
  expect_equal(length(bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE, check.chr.name=FALSE)), 2)
  bin.df$chr = 180
  suppressWarnings(expect_error(bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE), "Couldn't guess"))
  suppressWarnings(expect_error(bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE, no.checks=TRUE)))
  bin.df$chr = 18
  expect_equal(length(bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE)), 2)
  suppressWarnings(expect_error(bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE, check.chr.name=FALSE)))
})

test_that("Works with different chunk sizes",{
  expect_equal(length(bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE, chunk.size=1e3)), 2)
  expect_equal(length(bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE, chunk.size=1e5)), 2)
  expect_equal(length(bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE, chunk.size=1e7)), 2)
})

test_that("Counts concordant reads",{
  bb = bin.bam(bam.f, bin.df, appendIndex.outfile=FALSE)
  expect_equal(any(bb$bc$bc!=0), TRUE)
  expect_equal(all(tail(bb$bc$bc)==0), TRUE)
})

test_that("Writes files",{
  bb = bin.bam(bam.f, bin.df, outfile.prefix = "test.tsv", appendIndex.outfile=FALSE)
  expect_equal(file.exists("test.tsv"), TRUE)
  file.remove("test.tsv")
  bb = bin.bam(bam.f, bin.df, outfile.prefix = "test.tsv", appendIndex.outfile=TRUE)
  expect_equal(file.exists(paste0("test.tsv",c(".bgz",".bgz.tbi"))),c(TRUE,TRUE))
  file.remove(paste0("test.tsv",c(".bgz",".bgz.tbi")))
})

test_that("Stops if no reads found or not",{
  expect_error(bin.bam(bam.f, tail(bin.df), appendIndex.outfile=FALSE), "no reads found")
  expect_equal(length(bin.bam(bam.f, tail(bin.df), appendIndex.outfile=FALSE, no.checks=TRUE)), 2)
})
