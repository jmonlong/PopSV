##' Find random genomic regions with the same properties as a set of input regions. The properties are defined using annotation tracks. For example, we can simulate random regions that overlap genes as much as the input regions. It can be used to perform enrichment analysis while controling for specific properties.
##' @title Select control genomic regions for enrichment analysis
##' @param cnv.gr the input regions
##' @param feat.grl a list of the GRanges defining the features to fit. E.g. gene annotation, centromere, etc
##' @param nb.class the number of size class to speed up the computation. Default is 20.
##' @param nb.cores the number of cores to use. Default is 3.
##' @param redo.duplicates should duplicate regions be reselected. Default is TRUE.
##' @param seed.nb.max the maximum number of genomic position to use as seeds. Default is 1e5.
##' @param min.nb.gr the minimum number of control regions. If NULL (default), as many controls as input are simulated.
##' @param chr.prefix the chromosome name prefix. Default is "" (no prefix). Other value could be "chr" if chromosome are defined as 'chr1', 'chr2', etc.
##' @return a GRanges object defining the control regions
##' @author Jean Monlong
##' @import GenomicRanges
##' @import GenomeInfoDb
##' @export
draw.controls <- function(cnv.gr, feat.grl, nb.class=20, nb.cores=3, redo.duplicates=TRUE, seed.nb.max=1e5, min.nb.gr=NULL, chr.prefix=""){

  randGR.bp <- function(n=10){
    seql.1.22 = seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[paste0("chr",1:22)]
    chrs = apply(rmultinom(n, 1, (seql.1.22/1e3)/sum(seql.1.22/1e3)),2,function(e)which(e==1))
    starts = runif(n, 0, seql.1.22[chrs])
    return(GenomicRanges::GRanges(paste0(chr.prefix, chrs), IRanges::IRanges(starts, width=1)))
  }

  if(is.data.frame(cnv.gr)){
    cnv.gr = GenomicRanges::makeGRangesFromDataFrame(cnv.gr)
  }
  if(!is.null(min.nb.gr) && min.nb.gr > length(cnv.gr)){
    cnv.gr = cnv.gr[sample.int(length(cnv.gr), min.nb.gr, replace=TRUE)]
  }
  if(is.list(feat.grl) && is.data.frame(feat.grl[[1]])){
    feat.grl = lapply(feat.grl, GenomicRanges::makeGRangesFromDataFrame)
  }
  d.gr = randGR.bp(max(seed.nb.max, 3*length(cnv.gr)))
  d.l = parallel::mclapply(feat.grl, function(feat.gr){
    dtn = GenomicRanges::distanceToNearest(d.gr, feat.gr)
    as.data.frame(dtn)$distance
  }, mc.cores=nb.cores)
  d.df = as.data.frame(d.l)
  ol.l = parallel::mclapply(feat.grl, function(feat.gr)GenomicRanges::overlapsAny(cnv.gr, feat.gr), mc.cores=nb.cores)
  if(is.null(cnv.gr$sample)) {
    cnv.gr$sample = ""
  }
  cnv.df = data.frame(as.data.frame(ol.l), sample=cnv.gr$sample)
  ol.prof = unique(cnv.df[,1:length(feat.grl), drop=FALSE])
  null.gr = parallel::mclapply(1:nrow(ol.prof), function(ii){
    gr.ii = cnv.gr[which(apply(t(as.matrix(cnv.df[,1:length(feat.grl), drop=FALSE]))==unlist(ol.prof[ii,]),2,all))]
    widths = width(gr.ii)
    w.class = cut(widths, breaks=c(0,unique(quantile(widths, probs=ppoints(nb.class))), Inf))
    w.class = factor(w.class, levels=unique(w.class))
    tt = tapply(1:length(gr.ii), w.class, function(iii){
      w.i = widths[iii]
      demi.w = mean(w.i, na.rm=TRUE) / 2
      good.d = apply(t(as.matrix(d.df)-demi.w)*ifelse(unlist(ol.prof[ii,]),1,-1),2, function(x)all(x<0))
      gr = GenomicRanges::resize(d.gr[sample(which(good.d), length(w.i), TRUE)], w.i, fix="center")
      gr$sample = gr.ii$sample[iii]
      gr
    })
    res = unlist(do.call(GenomicRanges::GRangesList, tt[which(!unlist(lapply(tt,is.null)))]))
    names(res) = NULL
    res
  }, mc.cores=nb.cores)
  null.gr = unlist(do.call(GenomicRanges::GRangesList, null.gr))
  if(all(unique(null.gr$sample)=="")) {
    null.gr$sample = NULL
  }
  if(redo.duplicates && any((dup = duplicated(as.data.frame(null.gr))))){
    redo.gr = draw.controls(null.gr[which(dup)], feat.grl, nb.class=nb.class, nb.cores=nb.cores, redo.duplicates=FALSE, min.nb.gr=NULL, chr.prefix = chr.prefix)
    null.gr = c(null.gr[which(!dup)], redo.gr)
  }
  null.gr
}
