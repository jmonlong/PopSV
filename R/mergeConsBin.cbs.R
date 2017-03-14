##' Circular Binary Segmentation of the Z-scores. (NOT TESTED WELL YET)
##' @title Segmentation using CBS (NOT READY)
##' @param df a data.frame with at least 'chr', 'start' and 'end' columns.
##' @param pv.th the P-value threshold used for the segmentation and the segment filtering. Default is 0.01.
##' @param stitch.dist the stitching distance, i.e. the maximum distance at which two bins will be merged. 
##' @return a data.frame, similar to the input but with merged rows.
##' @author Jean Monlong
##' @keywords internal
mergeConsBin.cbs <- function(df, pv.th=.01, stitch.dist = 10000) {
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

  cna.l = lapply(unique(df$chr), function(chr.i){
    chr.ii = which(df$chr==chr.i)
    cna.o = DNAcopy::CNA(sign(df$z[chr.ii])*log10(df$qv[chr.ii]), df$chr[chr.ii], df$start[chr.ii])
    cna.s = DNAcopy::segment(cna.o, alpha=pv.th, verbose=0)
    cna.s$output
  })
  cna.s = as.data.frame(data.table::rbindlist(cna.l))
  cna.s = cna.s[which(abs(cna.s$seg.mean) > abs(log10(pv.th))),]
  
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
