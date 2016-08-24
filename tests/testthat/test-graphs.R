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

