##' Merge nearby bins using a user-defined "stitching" distance. Duplication and deletions are stitched separately.
##' @title Merge nearby bins
##' @param df a data.frame with at least 'chr', 'start' and 'end' columns.
##' @param col.mean the columns from which the average across the bins to merge will be computed.
##' @param stitch.dist the stitching distance, i.e. the maximum distance at which two bins will be merged. 
##' @return a data.frame, similar to the input but with merged rows. 
##' @author Jean Monlong
##' @export
mergeConsBin.reduce <- function(df, col.mean=c("z","pv","qv","fc"), stitch.dist=1e4){
  if(nrow(df)==0) return(df)
  
  fun3 <- function(x, FUN=mean){
    if(length(x)>2){
      FUN(x[c(2,length(x))])
    } else if(length(x)==2){
      return(FUN(x))
    } else {
      return(x)
    }
  }
  
  gr.f = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), z=z))
  dup.gr = with(df[which(df$z>0),], GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
  del.gr = with(df[which(df$z<0),], GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
  df$red.i = NA

  ## Merge duplications 
  gr.red.dup = GenomicRanges::reduce(dup.gr, min.gapwidth=stitch.dist)
  ol.dup = GenomicRanges::findOverlaps(gr.f, gr.red.dup)
  df$red.i[IRanges::queryHits(ol.dup)] = paste0("dup",IRanges::subjectHits(ol.dup))
  ## Merge deletions
  gr.red.del = GenomicRanges::reduce(del.gr, min.gapwidth=stitch.dist)
  ol.del = GenomicRanges::findOverlaps(gr.f, gr.red.del)
  df$red.i[IRanges::queryHits(ol.del)] = paste0("del",IRanges::subjectHits(ol.del))

  merge.event.f <- function(df.f){
    df.o = with(df.f, data.frame(start=min(start), end=max(end), nb.bin.cons=nrow(df.f)))
    cbind(df.o,t(apply(df.f[,intersect(colnames(df.f),col.mean),drop=FALSE],2,fun3)))
  }

  red.i = chr = . = NULL ## Uglily appease R checks
  df.o = as.data.frame(dplyr::do(dplyr::group_by(df, red.i, chr), merge.event.f(.)))
  df.o$red.i = NULL
  df.o
}
