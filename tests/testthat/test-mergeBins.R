context("Merging consecutive bins with abnormal coverage")

## Creating a data.frame with bins, z scores and abnormal regions
z.df = fragment.genome.hg19(5e5)
z.df$z = rnorm(nrow(z.df))
for(ii in 1:10){
  cnv.s = round(runif(1,1,nrow(z.df)))
  cnv.w = round(runif(1,1,min(c(10,nrow(z.df)-cnv.s))))
  z.df$z[cnv.s:(cnv.s+cnv.w-1)] = sample(c(-1,1),1)*runif(cnv.w, -5,100)
}
## Dummy information
z.df$fc = runif(nrow(z.df),0,5)
z.df$mean.cov = rnorm(nrow(z.df),3000,500)
## P-values
fdr = fdrtool.quantile(z.df$z,quant.int=seq(.6,.99,.01))
z.df$qv = z.df$pv = NA
z.df$pv[which(z.df$z > 0)] = 2 * stats::pnorm(-abs(z.df$z[which(z.df$z > 0)]), 0, fdr$sigma.est.dup)
z.df$pv[which(z.df$z < 0)] = 2 * stats::pnorm(-abs(z.df$z[which(z.df$z < 0)]), 0, fdr$sigma.est.del)
if (any(z.df$pv == 0, na.rm = TRUE)){
  z.df$pv[which(z.df$pv == 0)] = .Machine$double.xmin
}
z.df$qv = stats::p.adjust(z.df$pv, method = "fdr")


test_that("CBS output looks good", {
  z.m = mergeConsBin.cbs(z.df)
  expect_true(nrow(z.m)>0)
  expect_true(length(unique(z.m$end-z.m$start))>1)
  expect_true(all(c("pv","qv","mean.cov","fc") %in% colnames(z.m)))
})

test_that("Reduce merge output looks good", {
  z.m = mergeConsBin.reduce(subset(z.df, qv<.01))
  expect_true(nrow(z.m)>0)
  expect_true(length(unique(z.m$end-z.m$start))>1)
  expect_true(sum((z.m$end-z.m$start)/1e3)/sum((z.df$end-z.df$start)/1e3)<.1)
  expect_true(all(c("pv","qv","mean.cov","fc") %in% colnames(z.m)))
})

test_that("Z-pairs merge output looks good", {
  z.m = mergeConsBin.z(z.df, fdr.th = .01, sd.null = 1, nb.sim=1e3)
  expect_true(nrow(z.m)>0)
  expect_true(length(unique(z.m$end-z.m$start))>1)
  expect_true(sum((z.m$end-z.m$start)/1e3)/sum((z.df$end-z.df$start)/1e3)<.1)
  expect_true(all(c("pv","qv","mean.cov","fc") %in% colnames(z.m)))
})

test_that("Fragmented calls recovered", {
  z.df$z = rnorm(nrow(z.df))
  cnv.s = round(runif(1,1,nrow(z.df)-20))
  while(length(unique(z.df$chr[(cnv.s-5):(cnv.s+20)]))>1){
    cnv.s = round(runif(1,1,nrow(z.df)-20))
  }
  cnv.w1 = round(runif(1,1,10))
  cnv.w2 = round(runif(1,1,10))
  z.df$z[cnv.s:(cnv.s+cnv.w1-1)] = runif(cnv.w1, 10,100)
  z.df$z[(cnv.s+cnv.w1+1):(cnv.s+cnv.w1+cnv.w2)] = runif(cnv.w2, 10,100)
  fdr = fdrtool.quantile(z.df$z,quant.int=seq(.6,.99,.01))
  z.df$qv = z.df$pv = NA
  z.df$pv[which(z.df$z > 0)] = 2 * stats::pnorm(-abs(z.df$z[which(z.df$z > 0)]), 0, fdr$sigma.est.dup)
  z.df$pv[which(z.df$z < 0)] = 2 * stats::pnorm(-abs(z.df$z[which(z.df$z < 0)]), 0, fdr$sigma.est.del)
  if (any(z.df$pv == 0, na.rm = TRUE)){
    z.df$pv[which(z.df$pv == 0)] = .Machine$double.xmin
  }
  z.df$qv = stats::p.adjust(z.df$pv, method = "fdr")
  z.m = mergeConsBin.reduce(subset(z.df, qv<.01))
  expect_true(sum(z.m$chr==z.df$chr[cnv.s] & z.m$start>z.df$start[cnv.s-2] & z.m$start<z.df$start[cnv.s+cnv.w1+cnv.w2+1])>1)
  z.m = mergeConsBin.reduce(subset(z.df, qv<.01), stitch.dist = 1e6)
  expect_true(sum(z.m$chr==z.df$chr[cnv.s] & z.m$start>z.df$start[cnv.s-2] & z.m$start<z.df$start[cnv.s+cnv.w1+cnv.w2+1])==1)
  z.m = mergeConsBin.z(z.df, stitch.dist = 1e6, nb.sim=1e3)
  expect_true(sum(z.m$chr==z.df$chr[cnv.s] & z.m$start>z.df$start[cnv.s-2] & z.m$start<z.df$start[cnv.s+cnv.w1+cnv.w2+1])==1)
})
