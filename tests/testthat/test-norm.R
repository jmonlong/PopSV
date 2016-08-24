context("Test RD normalization")

## Create fake bins and bin counts
bins.df = data.frame(chr=sample.int(22, 1000, TRUE), start=seq(1,1e4,10))
bins.df$end = bins.df$start + 9
bc.mat = replicate(10, abs(rnorm(nrow(bins.df), 100, 10)))
colnames(bc.mat) = paste0("s",1:10)
bc.df = cbind(bins.df, bc.mat)

ns.f = "tempmsd.tsv"

test_that("Targeted normalization works, and with different parameters", {
              bins.df = chunk.bin(bins.df, bg.chunk.size=200, sm.chunk.size=20)
              expect_true(all(c("bg.chunk","sm.chunk","bin") %in% colnames(bins.df)))
              bins.sub = subset(bins.df, sm.chunk==bins.df$sm.chunk[1])
              tn.o = tn.norm(bc.df, "s1", nb.support.bins = 20, bins=bins.sub$bin)
              expect_true(any(!is.na(tn.o$norm.stats[,-(1:3)])))
              expect_true(any(!is.na(tn.o$bc.norm[,-(1:3)])))
              tn.o = tn.norm(bc.df, "s1", nb.support.bins = 20, bins=bins.sub$bin, force.diff.chr = FALSE)
              expect_true(any(!is.na(tn.o$norm.stats[,-(1:3)])))
              expect_true(any(!is.na(tn.o$bc.norm[,-(1:3)])))
              tn.o = tn.norm(bc.df, "s1", nb.support.bins = 20, bins=bins.sub$bin, norm="trim")
              expect_true(any(!is.na(tn.o$norm.stats[,-(1:3)])))
              expect_true(any(!is.na(tn.o$bc.norm[,-(1:3)])))
              ns.df = tn.o$norm.stats
              write.table(ns.df, file=ns.f, sep="\t", row.names=FALSE, quote=FALSE)
          })

test_that("Targeted normalization QC works", {
              graphs = tn.norm.qc(ns.f, out.pdf="temp.pdf")
              expect_true(file.remove("temp.pdf"))
              graphs = tn.norm.qc(ns.f, out.pdf="temp.pdf", bin.size=TRUE)
              expect_true(file.remove("temp.pdf"))
              graphs = tn.norm.qc.div(ns.f, out.pdf="temp.pdf", chunk.size=30)
              expect_true(file.remove("temp.pdf"))
              expect_true(file.remove(ns.f))
          })

test_that("Other normalization works", {
              norm.o = med.norm(bc.df)
              expect_true(any(!is.na(norm.o$norm.stats[,-(1:3)])))
              expect_true(any(!is.na(norm.o$bc.norm[,-(1:3)])))
              norm.o = medvar.norm(bc.df, colnames(bc.mat))
              expect_true(any(!is.na(norm.o$norm.stats[,-(1:3)])))
              expect_true(any(!is.na(norm.o$bc.norm[,-(1:3)])))
              norm.o = medvar.norm(bc.df, colnames(bc.mat), z.poisson=TRUE)
              expect_true(any(!is.na(norm.o$norm.stats[,-(1:3)])))
              expect_true(any(!is.na(norm.o$bc.norm[,-(1:3)])))
              norm.o = quant.norm(bc.df)
              expect_true(any(!is.na(norm.o$norm.stats[,-(1:3)])))
              expect_true(any(!is.na(norm.o$bc.norm[,-(1:3)])))
              norm.o = pca.norm(bc.df)
              expect_true(any(!is.na(norm.o$norm.stats[,-(1:3)])))
              expect_true(any(!is.na(norm.o$bc.norm[,-(1:3)])))
              norm.o = tnK.norm(bc.df, "s1", cont.sample="s2")
              expect_true(any(!is.na(norm.o[,-(1:3)])))
          })
