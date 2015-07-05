##' Circular Binary Segmentation of the Z-scores.
##' @title Segmentation using CBS
##' @param df a data.frame with at least 'chr', 'start' and 'end' columns.
##' @param pv.th the P-value threshold used for the segmentation and the segment filtering. Default is 0.01.
##' @return a data.frame, similar to the input but with merged rows.
##' @author Jean Monlong
##' @export
mergeConsBin.cbs <- function(df, pv.th=.01) {
  if (nrow(df) == 0)
    return(df)

  col.mean = c("z", "fc", "mean.cov", "GCcontent", "lowC", "map")
  col.mean.log = c("pv", "qv")

  fun3 <- function(x, FUN = mean, log.x=FALSE) {
    if(log.x){x = log(x)}
    if (length(x) > 2) {
      res.x = FUN(x[2:(length(x)-1)])
    } else if (length(x) == 2) {
      res.x = FUN(x)
    } else {
      res.x = x
    }
    if(log.x){res.x = exp(res.x)}
    res.x
  }

  cna.o = with(df, DNAcopy::CNA(z, chr, start))
  cna.s = DNAcopy::segment(cna.o, alpha=pv.th, undo.splits="prune")

  gr.f = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), z = z))
  dup.gr = with(df[which(df$z > 0), ], GenomicRanges::GRanges(chr, IRanges::IRanges(start,
    end)))
  del.gr = with(df[which(df$z < 0), ], GenomicRanges::GRanges(chr, IRanges::IRanges(start,
    end)))
  df$red.i = NA

  ## Merge duplications
  gr.seg = with(cna.s$output, GenomicRanges::GRanges(chrom, IRanges::IRanges(loc.start, loc.end)))
  ol.o = GenomicRanges::findOverlaps(gr.f, gr.seg)
  df$red.i[IRanges::queryHits(ol.o)] = paste0("seg", IRanges::subjectHits(ol.o))

  merge.event.f <- function(df.f) {
    df.o = with(df.f, data.frame(start = min(start), end = max(end), nb.bin.cons = nrow(df.f)))
    cbind(df.o,
          t(apply(df.f[, intersect(colnames(df.f), col.mean), drop = FALSE], 2, fun3)),
          t(apply(df.f[, intersect(colnames(df.f), col.mean.log), drop = FALSE], 2, fun3, log.x=TRUE))
          )
  }

  red.i = chr = . = NULL  ## Uglily appease R checks
  df.o = as.data.frame(dplyr::do(dplyr::group_by(df, red.i, chr), merge.event.f(.)))
  df.o$red.i = NULL
  df.o[which(df.o$pv<pv.th),]
}
