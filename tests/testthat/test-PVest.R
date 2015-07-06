context("P-value estimation from the Z-score distribution")

## Creating a data.frame with bins, z scores and abnormal regions
z.df = data.frame(chr=sample(1:2,1e4, replace=TRUE), start=sample(seq(1,1e7,1e4),1e4,TRUE), stringsAsFactors=FALSE)
z.df$end = z.df$start + 1e4 -1
z.df$z = rnorm(nrow(z.df))

flatness <- function(x){
  x.t = table(cut(x,10))
  (max(x.t)-min(x.t))/length(x)
}

test_that("Normal Z-scores leads to flat P-values using quantile approach",{
  fdr = fdrtool.quantile(z.df$z,quant.int=seq(.9,.99,.01), plot = FALSE)
  expect_true(flatness(fdr$pval)<.1)
  fdr = fdrtool.quantile(z.df$z,quant.int=seq(.6,.99,.01), plot = FALSE)
  expect_true(flatness(fdr$pval)<.1)
  fdr = fdrtool.quantile(z.df$z,quant.int=seq(.3,.99,.01), plot = FALSE)
  expect_true(flatness(fdr$pval)<.1)
})

test_that("Normal Z-scores + outliers leads to flat P-values using quantile approach",{
  z.df$z[sample.int(nrow(z.df),100)] = runif(100,-100,100)
  fdr = fdrtool.quantile(z.df$z,quant.int=seq(.9,.99,.01), plot = FALSE)
  expect_true(flatness(fdr$pval)<.1)
  fdr = fdrtool.quantile(z.df$z,quant.int=seq(.6,.99,.01), plot = FALSE)
  expect_true(flatness(fdr$pval)<.1)
  fdr = fdrtool.quantile(z.df$z,quant.int=seq(.3,.99,.01), plot = FALSE)
  expect_true(flatness(fdr$pval)<.1)
})

test_that("Normal Z-scores leads to flat P-values using quantile 2N approach",{
  fdr = fdrtool.quantile.2N(z.df$z,plot = FALSE, min.prop.null = .9)
  expect_true(flatness(fdr$pval)<.1)
  fdr = fdrtool.quantile.2N(z.df$z, plot = FALSE, min.prop.null = .6)
  expect_true(flatness(fdr$pval)<.1)
  fdr = fdrtool.quantile.2N(z.df$z, plot = FALSE, min.prop.null = .4)
  expect_true(flatness(fdr$pval)<.1)
})

test_that("Normal Z-scores + outliers leads to flat P-values using quantile 2N approach",{
  z.df$z[sample.int(nrow(z.df),100)] = runif(100,-100,100)
  fdr = fdrtool.quantile.2N(z.df$z,plot = FALSE, min.prop.null = .9)
  expect_true(flatness(fdr$pval)<.1)
  fdr = fdrtool.quantile.2N(z.df$z, plot = FALSE, min.prop.null = .6)
  expect_true(flatness(fdr$pval)<.1)
  fdr = fdrtool.quantile.2N(z.df$z, plot = FALSE, min.prop.null = .4)
  expect_true(flatness(fdr$pval)<.1)
})

