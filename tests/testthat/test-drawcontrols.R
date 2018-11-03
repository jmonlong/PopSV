context("Draw controls")

## Create inputs
cnv.df = data.frame(chr=1, start=round(runif(10,0,1e3)))
cnv.df$end = cnv.df$start + runif(10,10,100)
cnv.gr = with(cnv.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start=start, end=end)))

feat.grl = list(feat1=GenomicRanges::GRanges(1, IRanges::IRanges(start=round(runif(100,0,1e3)), width=runif(100,10,20))), feat2=GenomicRanges::GRanges(1, IRanges::IRanges(start=round(runif(10,0,1e3)), width=runif(10,10,20))))

d.gr = GenomicRanges::GRanges(1, IRanges::IRanges(start=round(runif(10,0,1e3)), width=10))

test_that("It runs with different parameters",{
              c.gr = draw.controls(cnv.gr, feat.grl, nb.cores=1)
              expect_true(all(sort(GenomicRanges::width(c.gr)) == sort(GenomicRanges::width(cnv.gr))))
              expect_equal(length(c.gr), length(cnv.gr))
              c.gr = draw.controls(cnv.df, feat.grl, nb.cores=1, nb.class=2, redo.duplicates = FALSE)
              expect_true(all(sort(GenomicRanges::width(c.gr)) == sort(GenomicRanges::width(cnv.gr))))
              expect_equal(length(c.gr), length(cnv.gr))
              c.gr = draw.controls(cnv.gr, feat.grl, nb.cores=1, dist.gr=d.gr, redo.duplicates = FALSE)
              expect_true(all(sort(GenomicRanges::width(c.gr)) == sort(GenomicRanges::width(cnv.gr))))
              expect_equal(length(c.gr), length(cnv.gr))
          })

