context("Test cases: RD normalization + Z-scores")

## Create fake bins and bin counts
bins.df = data.frame(chr=sample.int(22, 1000, TRUE), start=seq(1,1e4,10))
bins.df$end = bins.df$start + 9
bc.mat = replicate(10, rnorm(nrow(bins.df), 100, 10))
colnames(bc.mat) = paste0("s",1:10)

## Write files
files.df = data.frame(sample=colnames(bc.mat), stringsAsFactors=FALSE)
files.df$bc.gc.gz = paste0(files.df$sample, "-test-bc.tsv")
files.df$bc.gc.norm = paste0(files.df$sample, "-test-bc-norm.tsv")
files.df$z = paste0(files.df$sample, "-test-z.tsv")
files.df$fc = paste0(files.df$sample, "-test-fc.tsv")
dump = sapply(1:ncol(bc.mat), function(ii){
                  bins.df$bc = bc.mat[,ii]
                  write.table(bins.df, file=files.df$bc.gc.gz[ii], sep="\t", row.names=FALSE, col.names=TRUE)
              })

bc.f = "tempbc.tsv"
bc.df = cbind(bins.df, bc.mat)
write.table(bc.df, file=bc.f, sep="\t", row.names=FALSE, col.names=TRUE)

ns.f = "tempns.tsv"

exp.outfiles = c(paste0(c(files.df$z[2], files.df$fc[2], files.df$bc.gc.norm[2]), ".bgz"),
paste0(c(files.df$z[2], files.df$fc[2], files.df$bc.gc.norm[2]), ".bgz.tbi"))

test_that("Normalization works, and with different parameters", {
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

test_that("Case z-score runs with different parameters",{
              res.df = tn.test.sample("s2", files.df, "s1", norm.stats.f = ns.f)
              expect_true(all(file.remove(exp.outfiles)))
              res.df = tn.test.sample("s2", files.df, "s1", bc.ref.f=bc.f, z.poisson=TRUE, aberrant.cases = TRUE,  norm.stats.f = ns.f)
              expect_true(all(file.remove(exp.outfiles)))
          })


## Remove files
dump = file.remove(files.df$bc.gc.gz, bc.f, ns.f)



