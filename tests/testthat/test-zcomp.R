context("Z-score computation")

## Create fake bins and bin counts
bc.df = data.frame(chr=1, start=seq(1,1e4,10))
bc.df$end = bc.df$start + 9
bc.mat = replicate(10, rnorm(nrow(bc.df), 100, 10))
colnames(bc.mat) = paste0("s",1:10)

## Write files
files.df = data.frame(sample=colnames(bc.mat), stringsAsFactors=FALSE)
files.df$bc.gc.norm.gz = paste0(files.df$sample, "-test-bc.tsv")
files.df$z = paste0(files.df$sample, "-test-z.tsv")
files.df$fc = paste0(files.df$sample, "-test-fc.tsv")
dump = sapply(1:10, function(ii){
                  bc.df$bc = bc.mat[,ii]
                  write.table(bc.df, file=files.df$bc.gc.norm.gz[ii], sep="\t", row.names=FALSE, col.names=TRUE)
              })

bc.f = "tempbc.tsv"
bc.df = cbind(bc.df, bc.mat)
write.table(bc.df, file=bc.f, sep="\t", row.names=FALSE, col.names=TRUE)

test_that("It runs with different parameters",{
              z.df = z.comp(bc.f, files.df)
              expect_true(all(file.exists(paste0(files.df$z, ".bgz"))))
              expect_true(all(file.exists(paste0(files.df$fc, ".bgz"))))
              expect_true(file.exists("ref-msd.tsv"))
              z.df = z.comp(bc.f, files.df, z.poisson=TRUE)
              expect_true(all(file.exists(paste0(files.df$z, ".bgz"))))
              expect_true(all(file.exists(paste0(files.df$fc, ".bgz"))))
              expect_true(file.exists("ref-msd.tsv"))
              z.df = z.comp(bc.f, files.df, ref.samples=head(files.df$sample))
              expect_true(all(file.exists(paste0(files.df$z, ".bgz"))))
              expect_true(all(file.exists(paste0(files.df$fc, ".bgz"))))
              expect_true(file.exists("ref-msd.tsv"))
              z.df = z.comp(files.df=files.df)
              expect_true(all(file.exists(paste0(files.df$z, ".bgz"))))
              expect_true(all(file.exists(paste0(files.df$fc, ".bgz"))))
              expect_true(file.exists("ref-msd.tsv"))
              z.df = z.comp(bc.f, files.df=files.df, chunk.size=100)
          })

test_that("Bad input raise errors",{
              expect_error(z.comp("badpath", files.df), "file not found")
              expect_error(z.comp(files.df=files.df, files.col=NULL), "non-NULL")
              expect_error(z.comp(bc.f, files.df=files.df, ref.samples = c("badsample",files.df$sample)), "not present")
          })


## Remove files
dump = file.remove(files.df$bc.gc.norm.gz, paste0(files.df$z, ".bgz"), paste0(files.df$z, ".bgz.tbi"), paste0(files.df$fc, ".bgz"), paste0(files.df$fc, ".bgz.tbi"), bc.f, "ref-msd.tsv")
