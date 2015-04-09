context("Frequency computation from ranges")

bin.df = data.frame(chr=sample(1:2,1e2, replace=TRUE), start=sample(seq(1,1e6,1e4)), stringsAsFactors=FALSE)
bin.df$end = bin.df$start + sample(1:5, 1e2, replace=TRUE) * 1e4 -1

test_that("Missing columns tigger an error", {
  bad.df = bin.df[,c("chr","start")]
  expect_error(freq.range(bad.df), "Missing column")
  bad.df = bin.df[,c("end","start")]
  expect_error(freq.range(bad.df), "Missing column")
  bad.df = bin.df[,c("chr","end")]
  expect_error(freq.range(bad.df), "Missing column")
  good.df = bin.df[,c("chr","end", "start")]
  expect_true(nrow(freq.range(good.df))>0)
})

test_that("The total amount of caled region is conserved", {
  fr.df = freq.range(bin.df)
  expect_true(nrow(fr.df)>0)
  expect_equal(sum(bin.df$end-bin.df$start+1), sum((fr.df$end-fr.df$start+1)*fr.df$nb))
})

test_that("Graphs doesn't throw an error", {
  pdf("tmp.pdf")
  fr.df = freq.range(bin.df, plot=TRUE)
  expect_true(nrow(fr.df)>0)
  dev.off()
  file.remove("tmp.pdf")
})


