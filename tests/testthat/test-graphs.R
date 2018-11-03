context("Graphs")

## Create calls
res.df = data.frame(chr=sample.int(22, 100, TRUE), start=sample(seq(1,1e4,10), 100))
res.df$end = res.df$start + 9
res.df$z = rnorm(nrow(res.df))
res.df$pv = runif(nrow(res.df), 0, .01)
res.df$qv = res.df$pv*2
res.df$fc = rnorm(nrow(res.df), sample(c(0,1,3,4,5), nrow(res.df), replace=TRUE), .1)
res.df$nb.bin.cons = sample(1:10, nrow(res.df), replace=TRUE)
res.df$cn2.dev = runif(nrow(res.df))
res.df$sample = sample(paste0("samp",1:4), nrow(res.df), replace=TRUE)

pdf.f = "temp.pdf"

test_that("Summary graphs", {
              graphs = sv.summary(res.df, out.pdf=pdf.f, print=FALSE)
              expect_true(file.remove(pdf.f))
              expect_true(length(graphs)>0)
              expect_true(length(graphs$graphs.l)>0)
          })

test_that("Chromosome plots", {
              graphs = chrplot(res.df, bin.size=3e6)
              expect_true(length(graphs)>0)
              graphs = chrplot(res.df, bin.size=3e6, showSampleNames=TRUE)
              expect_true(length(graphs)>0)
          })

## Fake output files
bins.df = data.frame(chr=1, start=seq(1,1e4,10))
bins.df$end = bins.df$start + 9
bc.mat = replicate(10, abs(rnorm(nrow(bins.df), 100, 10)))
colnames(bc.mat) = paste0("s",1:10)
bc.df = cbind(bins.df, bc.mat)
bc.f = "tempbc.tsv"
write.table(bc.df, file=bc.f, row.names=FALSE, sep="\t", quote=FALSE)

msd.df = data.frame(bins.df[sample.int(nrow(bins.df)),], m=rnorm(nrow(bins.df), 100, 10), sd=rnorm(nrow(bins.df), 5, 1))
ns.f = "msd-test.tsv"
write.table(msd.df, file=ns.f, sep="\t", row.names=FALSE, col.names=TRUE)

res.df = rbind(res.df, data.frame(chr=1, start=110, end=130, z=10,pv=0, qv=0, fc=2, nb.bin.cons=2, cn2.dev=1, sample="s1"))

test_that("Example graph works with its many parameters", {
              comp.index.files(c(bc.f, ns.f), reorder=TRUE)
              pdf("temp.pdf")
              coverage.plot(1, 100, 130, bc.f, flanks=50)
              coverage.plot(1, 100, 130, bc.f, flanks=50, boxplot=TRUE)
              coverage.plot(1, 100, 130, bc.f, flanks=50, sv.df=res.df)
              coverage.plot(1, 100, 130, bc.f, flanks=50, sv.df=res.df, ref.samples=c("s2","s3","s4"))
              coverage.plot(1, 100, 130, bc.f, flanks=50, norm.stats.f=ns.f, sv.df=res.df)
              coverage.plot(1, 100, 130, bc.f, flanks=50, norm.stats.f=ns.f, samples="s1")
              coverage.plot(1, 100, 130, bc.f, flanks=50, norm.stats.f=ns.f, samples="s1", absolute.position =FALSE)
              dev.off()
              expect_true(file.remove("temp.pdf"))
          })

file.remove(paste0(bc.f, c(".bgz",".bgz.tbi")), paste0(ns.f, c(".bgz",".bgz.tbi")))
