context("Z-score computation")

## Create fake bins, bin counts and norm stats
bc.df = data.frame(chr=1, start=seq(1,1e4,10))
bc.df$end = bc.df$start + 9
bc.mat = replicate(10, rnorm(nrow(bc.df), 100, 10))
colnames(bc.mat) = paste0("s",1:10)
bc.mat = bc.mat[sample.int(nrow(bc.mat)),]

## Write files
files.df = data.frame(sample=colnames(bc.mat), stringsAsFactors=FALSE)
files.df$z = paste0(files.df$sample, "-test-z.tsv")
files.df$fc = paste0(files.df$sample, "-test-fc.tsv")

bc.f = "tempbc.tsv"
ns.f = "tempns.tsv"
bc.df = cbind(bc.df, bc.mat)
ns.df = bc.df[,1:3]
ns.df$m = rnorm(nrow(bc.df), 100, 10)
ns.df$sd = abs(rnorm(nrow(bc.df), 10, 10))
write.table(bc.df, file=bc.f, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(ns.df, file=ns.f, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

test_that("It runs with different parameters",{
              z.df = z.comp(bc.f, ns.f, files.df)
              expect_true(all(file.exists(paste0(files.df$z, ".bgz"))))
              expect_true(all(file.exists(paste0(files.df$fc, ".bgz"))))
              z.df = z.comp(bc.f, ns.f, files.df, z.poisson=TRUE)
              expect_true(all(file.exists(paste0(files.df$z, ".bgz"))))
              expect_true(all(file.exists(paste0(files.df$fc, ".bgz"))))
              z.df = z.comp(bc.f, ns.f, files.df, chunk.size=50)
              expect_true(all(file.exists(paste0(files.df$z, ".bgz"))))
              expect_true(all(file.exists(paste0(files.df$fc, ".bgz"))))
})

test_that("Bad input raise errors",{
  expect_error(z.comp("badpath", files.df), "file not found")
})


## Remove files
dump = file.remove(paste0(files.df$z, ".bgz"), paste0(files.df$z, ".bgz.tbi"), paste0(files.df$fc, ".bgz"), paste0(files.df$fc, ".bgz.tbi"), bc.f, ns.f)
