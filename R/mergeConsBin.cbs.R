##' Circular Binary Segmentation of the Z-scores. (NOT TESTED WELL YET)
##' @title Segmentation using CBS (NOT READY)
##' @param df a data.frame with at least 'chr', 'start' and 'end' columns.
##' @param pv.th the P-value threshold used for the segmentation and the segment filtering. Default is 0.01.
##' @param stitch.dist the stitching distance, i.e. the maximum distance at which two bins will be merged.
##' @param max.gap.size the maximum gap between bins allowed. Default is 100 kb. Calls will not span gap larger than this (e.g. centromere).
##' @return a data.frame, similar to the input but with merged rows.
##' @author Jean Monlong
##' @keywords internal
mergeConsBin.cbs <- function(df, pv.th=.01, stitch.dist = 10000, max.gap.size=1e5) {
  if (nrow(df) == 0)
    return(df)

  col.mean = c("z", "fc", "mean.cov", "GCcontent", "lowC", "map")
  col.min = c("pv", "qv")

  fun3 <- function(x, FUN3 = mean) {
    if (length(x) > 2) {
      res.x = FUN3(x[2:(length(x)-1)])
    } else if (length(x) == 2) {
      res.x = FUN3(x)
    } else {
      res.x = x
    }
    res.x
  }

  ## CBS
  if(!is.null(max.gap.size)){
    ## Create fake chromosome when there is a gap in the binned regions (larger than stitch.dist)
    seg = df$chr[-1] == df$chr[-nrow(df)] & df$start[-1] <= max.gap.size + df$end[-nrow(df)]
    chrs = rep(1:(1+sum(!seg)), diff(c(0, which(!seg), length(seg)+1)))
    cchrs = unique(cbind(chrs, df$chr))
    chrsTrans = cchrs[,2]
  } else {
    chrs = df$chr
  }
  cna.o = DNAcopy::CNA(sign(df$z)*log10(df$qv), chrs, df$start)
  cna.s = DNAcopy::segment(cna.o, alpha=pv.th, verbose=0)
  cna.s = as.data.frame(cna.s$output)
  cna.s = cna.s[which(abs(cna.s$seg.mean) > abs(log10(pv.th))),]
  if(!is.null(max.gap.size)){
    ## Back to real chromosomes
    cna.s$chrom = chrsTrans[cna.s$chrom]
  }

  ## Merge segments
  gr.f = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
  df$red.i = NA
  ## Duplications
  if(any(cna.s$seg.mean<0)){
    gr.seg = with(cna.s[which(cna.s$seg.mean<0),], GenomicRanges::GRanges(chrom, IRanges::IRanges(loc.start, loc.end)))
  } else {
    gr.seg = GenomicRanges::GRanges()
  }
  if(any(df$qv <= pv.th & df$z>0)){
    gr.seg = with(df[which(df$qv <= pv.th & df$z>0),], c(gr.seg, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))))
  }
  gr.seg = GenomicRanges::reduce(gr.seg, min.gapwidth = stitch.dist)
  ol.o = GenomicRanges::findOverlaps(gr.f, gr.seg)
  df$red.i[S4Vectors::queryHits(ol.o)] = paste0("dup", S4Vectors::subjectHits(ol.o))
  ## Deletions
  if(any(cna.s$seg.mean>0)){
    gr.seg = with(cna.s[which(cna.s$seg.mean>0),], GenomicRanges::GRanges(chrom, IRanges::IRanges(loc.start, loc.end)))
  } else {
    gr.seg = GenomicRanges::GRanges()
  }
  if(any(df$qv <= pv.th & df$z<0)){
    gr.seg = with(df[which(df$qv <= pv.th & df$z<0),], c(gr.seg, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))))
  }
  gr.seg = GenomicRanges::reduce(gr.seg, min.gapwidth = stitch.dist)
  ol.o = GenomicRanges::findOverlaps(gr.f, gr.seg)
  df$red.i[S4Vectors::queryHits(ol.o)] = paste0("del", S4Vectors::subjectHits(ol.o))

  merge.event.f <- function(df.f) {
    df.o = with(df.f, data.frame(start = min(start), end = max(end), nb.bin.cons = nrow(df.f)))
    cbind(df.o,
          t(apply(df.f[, intersect(colnames(df.f), col.mean), drop = FALSE], 2, fun3)),
          t(apply(df.f[, intersect(colnames(df.f), col.min), drop = FALSE], 2, fun3, FUN3=min))
          )
  }

  red.i = chr = . = NULL  ## Uglily appease R checks
  df = df[which(!is.na(df$red.i)),]
  df.o = as.data.frame(dplyr::do(dplyr::group_by(df, red.i, chr), merge.event.f(.)))
  df.o$red.i = NULL
  df.o
}
