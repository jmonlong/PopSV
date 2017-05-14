context("VCF convertion")

## Create inputs
cnv.df = data.frame(chr=1, start=round(runif(100, 0, 1e3)))
cnv.df$end = cnv.df$start + round(runif(100, 10, 100))
cnv.df$fc = runif(100, 0, 3)
cnv.df$nb.bin.cons = 1
cnv.df$sample = sample(letters[1:5], 100, TRUE)
vcf.temp = 'temp.vcf'

test_that("File is created and not empty", {
            expect_equal(writeVcf(cnv.df, vcf.temp), vcf.temp)
            vcf = read.table(vcf.temp)
            expect_true(nrow(vcf)>0)
            expect_true(file.remove(vcf.temp))
          })

test_that("Rerror is raised if wrong inputs", {
            cnv.df$nb.bin.cons = NULL
            expect_error(writeVcf(cnv.df, vcf.temp), "column")
          })

