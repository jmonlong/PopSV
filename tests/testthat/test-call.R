context("Abnormal coverage calling")

## Create fake bins, z, fc
bins.df = data.frame(chr=1, start=seq(1,1e4,10))
bins.df = rbind(bins.df, data.frame(chr=2, start=seq(1,1e4,10)))
bins.df$end = bins.df$start + 9
bins.df$GCcontent = runif(nrow(bins.df), .3,.7)
z.df = data.frame(bins.df, z=rnorm(nrow(bins.df)))
fc.df = data.frame(bins.df, fc=abs(rnorm(nrow(bins.df), 2, .01)))
msd.df = data.frame(bins.df[sample.int(nrow(bins.df)),], m=rnorm(nrow(bins.df), 1000, 10), sd=rnorm(nrow(bins.df), 10, 1))
ns.f = "msd-test.tsv"
write.table(msd.df, file=ns.f, sep="\t", row.names=FALSE, col.names=TRUE)

## Write files
samp.nocnv = "sampnocnv"
files.df = data.frame(sample=samp.nocnv, z=paste0(samp.nocnv, "-z.tsv"), fc=paste0(samp.nocnv, "-fc.tsv"), stringsAsFactors=FALSE)
write.table(z.df, file=files.df$z, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(fc.df, file=files.df$fc, sep="\t", row.names=FALSE, col.names=TRUE)

## Add some CNVs
pos1 = nrow(z.df)/3
pos2 = 2*pos1
z.df$z[pos1:(pos1+4)] = 10
z.df$z[pos2:(pos2+1)] = -10
fc.df$fc[pos1:(pos1+4)] = 4
fc.df$fc[pos2:(pos2+1)] = 0

## Write files
samp.cnv = "sampcnv"
files.df = rbind(files.df, data.frame(sample=samp.cnv, z=paste0(samp.cnv, "-z.tsv"), fc=paste0(samp.cnv, "-fc.tsv"), stringsAsFactors=FALSE))
write.table(z.df, file=files.df$z[2], sep="\t", row.names=FALSE, col.names=TRUE)
write.table(fc.df, file=files.df$fc[2], sep="\t", row.names=FALSE, col.names=TRUE)

test_that("It runs with different parameters",{
  res.df = call.abnormal.cov(files.df, samp.cnv, z.th="sdest", merge.cons.bins="stitch", out.pdf="test.pdf")
  res.df = call.abnormal.cov(files.df, samp.cnv, z.th="sdest", merge.cons.bins="stitch", out.pdf="test.pdf", stitch.dist=30)
  res.df = call.abnormal.cov(files.df, samp.cnv, z.th="sdest", merge.cons.bins="no", out.pdf="test.pdf")
  file.remove("test.pdf")
})

test_that("It runs with more different parameters",{
  res.df = call.abnormal.cov(files.df, samp.cnv, z.th="sdest", merge.cons.bins="zscores", out.pdf="test.pdf")
  res.df = call.abnormal.cov(files.df, samp.cnv, z.th="sdest", merge.cons.bins="cbs", out.pdf="test.pdf")
  file.remove("test.pdf")
})

test_that("It runs with even more different parameters",{
  res.df = call.abnormal.cov(files.df, samp.cnv, z.th="consbins", out.pdf="test.pdf")
  res.df = call.abnormal.cov(files.df, samp.cnv, z.th="sdest2N", out.pdf="test.pdf")
  file.remove("test.pdf")
})

test_that("It runs when there is no CNVs",{
  res.df = call.abnormal.cov(files.df, samp.nocnv, z.th="sdest", merge.cons.bins="stitch", out.pdf="test.pdf")
  res.df = call.abnormal.cov(files.df, samp.nocnv, z.th="sdest", merge.cons.bins="no", out.pdf="test.pdf")
  expect_true(file.remove("test.pdf"))
})

test_that("It still runs when there is no CNVs",{
              res.df = call.abnormal.cov(files.df, samp.nocnv, z.th="sdest", merge.cons.bins="zscores", out.pdf="test.pdf")
              expect_true(file.remove("test.pdf"))
  res.df = call.abnormal.cov(files.df, samp.nocnv, z.th="sdest", merge.cons.bins="cbs", out.pdf="test.pdf")
  expect_true(file.remove("test.pdf"))
})

test_that("It still still runs when there is no CNVs",{
  res.df = call.abnormal.cov(files.df, samp.nocnv, z.th="consbins", out.pdf="test.pdf")
  res.df = call.abnormal.cov(files.df, samp.nocnv, z.th="sdest2N", out.pdf="test.pdf")
  expect_true(file.remove("test.pdf"))
})


test_that("GC normalization runs",{
  res.df = call.abnormal.cov(files.df, samp.cnv, gc.df=bins.df)
})

test_that("Pvalues can be written",{
              res.df = call.abnormal.cov(files.df, samp.cnv, outfile.pv="temppv.tsv")
              pv.df = read.table("temppv.tsv.bgz", as.is=TRUE, header=TRUE)
              expect_true(nrow(pv.df)>0)
              expect_true(file.remove("test.pdf"))
})

test_that("Normalization stats are merged",{
  res.df = call.abnormal.cov(files.df, samp.cnv, norm.stats = ns.f)
  expect_true(any(colnames(res.df)=="mean.cov"))
})

test_that("Chrs can be excluded",{
  res.df = call.abnormal.cov(files.df, samp.cnv, aneu.chrs = 1)
  expect_true(all(res.df$chr!=1))
})

## Remove files
dump = file.remove(files.df$z, files.df$fc, ns.f, "temppv.tsv.bgz", "temppv.tsv.bgz.tbi")
