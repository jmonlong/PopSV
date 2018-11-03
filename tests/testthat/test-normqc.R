context("QC metrics on normalized bin count")

## Create fake bins and bin counts
bins.df = data.frame(chr=sample.int(22, 1000, TRUE), start=seq(1,1e4,10))
bins.df$end = bins.df$start + 9
bc.mat = replicate(10, rnorm(nrow(bins.df), 100, 10))
colnames(bc.mat) = paste0("s",1:10)
bc.df = cbind(bins.df, bc.mat)

test_that("It runs, and with different parameters", {
              res = normQC(bc.df, n.subset=50, win.size=10)
              nonna = lapply(res, function(x)any(!is.na(x)))
              expect_true(all(unlist(nonna)))
              pdf("temp.pdf")
              res = normQC(bc.df, n.subset=50, win.size=10, plot=TRUE)
              dev.off()
              expect_true(file.exists("temp.pdf"))
              file.remove("temp.pdf")
          })
