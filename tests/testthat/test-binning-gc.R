context("Automatic binning and GC content computation")


test_that("Different bin sizes works", {
  expect_true(nrow(fragment.genome.hp19(1e6))>0)
  expect_true(nrow(fragment.genome.hp19(1e5))>0)
  expect_true(nrow(fragment.genome.hp19(94862))>0)
})

test_that("Bins don't overlap", {
  gc.df = fragment.genome.hp19(94862)
  gc.gr = with(gc.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
  expect_equal(length(gc.gr), length(GenomicRanges::reduce(gc.gr, min.gapwidth=0)))
})

test_that("Only autosomal chromosomes are considered", {
  gc.df = fragment.genome.hp19(94862)
  expect_true(all(1:22 %in% unique(gc.df$chr)))
  expect_true(length(setdiff(unique(gc.df$chr), 1:22))==0)
})

test_that("Gets correct GC content", {
  ucsc.gc = read.table("../dataForTests/gc-hg19-chr18-30000-30100.bed.gz", skip=7, as.is=TRUE)
  colnames(ucsc.gc) = c("start", "GCcontent")
  ucsc.gc$chr = 18
  ucsc.gc$end = ucsc.gc$start + 4
  expect_equal(getGC.hg19(ucsc.gc[,c("chr","start","end")])$GCcontent , ucsc.gc$GCcontent/100, tolerance=.1)
})

test_that("Keeps bin order", {
  bin.df =  data.frame(chr=sample(1:22,1e2, replace=TRUE), start=seq(1,1e4,100), end=seq(99,9999,100))
  gc.df = getGC.hg19(bin.df)
  expect_equal(all(bin.df$chr==gc.df$chr & bin.df$start==gc.df$start), TRUE)
})

test_that("Checks for bins missing columns",{
  bin.df =  data.frame(chr=sample(1:22,1e2, replace=TRUE), start=seq(1,1e4,100), end=seq(99,9999,100))
  bad.df = bin.df[,c("chr","start")]
  expect_error( getGC.hg19(bad.df), "Missing column")
  bad.df = bin.df[,c("end","start")]
  expect_error( getGC.hg19(bad.df), "Missing column")
  bad.df = bin.df[,c("chr","end")]
  expect_error( getGC.hg19(bad.df), "Missing column")
  good.df = bin.df[,c("chr","end", "start")]
  expect_equal(nrow(getGC.hg19(good.df)), nrow(good.df))  
})

