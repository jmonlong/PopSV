##' Find random genomic regions with the same properties as a set of input regions. The properties are defined using annotation tracks. For example, we can simulate random regions that overlap genes as much as the input regions. It can be used to perform enrichment analysis while controling for specific properties.
##' @title Select control genomic regions for enrichment analysis
##' @param cnv.gr the input regions
##' @param feat.grl a list of the GRanges defining the features to fit. E.g. gene annotation, centromere, etc
##' @param nb.class the number of size class to speed up the computation. Default is 100. Inf will ensure exactly the same overlap proportion but might take some time if the size distribution is diverse.
##' @param nb.cores the number of cores to use. Default is 3.
##' @param redo.duplicates should duplicate regions be reselected. Default is TRUE.
##' @param seed.nb.max the maximum number of genomic position to use as seeds. Default is 1e5.
##' @param min.nb.gr the minimum number of control regions. If NULL (default), as many controls as input are simulated.
##' @param chr.prefix the chromosome name prefix. Default is "" (no prefix). Other value could be "chr" if chromosome are defined as 'chr1', 'chr2', etc.
##' @param dist.gr a GRanges defining the feature for which we want to control the distance to. Default is NULL, i.e. no control.
##' @return a GRanges object defining the control regions
##' @author Jean Monlong
##' @export
draw.controls <- function(cnv.gr, feat.grl, nb.class=100, nb.cores=3, redo.duplicates=TRUE, seed.nb.max=1e5, min.nb.gr=NULL, chr.prefix="", dist.gr=NULL){

  randGR.bp <- function(n=10){
    seql.1.22 = seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[paste0("chr",1:22)]
    chrs = apply(stats::rmultinom(n, 1, (seql.1.22/1e3)/sum(seql.1.22/1e3)),2,function(e)which(e==1))
    starts = stats::runif(n, 0, seql.1.22[chrs])
    return(GenomicRanges::GRanges(paste0(chr.prefix, chrs), IRanges::IRanges(starts, width=1)))
  }
  if(is.data.frame(cnv.gr)){
    cnv.gr = GenomicRanges::makeGRangesFromDataFrame(cnv.gr, keep.extra.columns = TRUE)
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
    ret = rep(stats::median(as.data.frame(dtn)$distance), length(d.gr))
    ret[S4Vectors::queryHits(dtn)] = as.data.frame(dtn)$distance
    ret
  }, mc.cores=nb.cores)
  d.df = as.data.frame(d.l)
  ol.l = parallel::mclapply(feat.grl, function(feat.gr)IRanges::overlapsAny(cnv.gr, feat.gr), mc.cores=nb.cores)
  if(is.null(cnv.gr$sample)) {
    cnv.gr$sample = ""
  }
  if(!is.null(dist.gr)){
    dtn = GenomicRanges::distanceToNearest(cnv.gr, dist.gr, ignore.strand=TRUE)
    cnv.gr$dist.feat = stats::median(as.data.frame(dtn)$distance)
    cnv.gr$dist.feat[S4Vectors::queryHits(dtn)] = as.data.frame(dtn)$distance
    dtn = GenomicRanges::distanceToNearest(d.gr, dist.gr)
    d.gr$dist.feat = stats::median(as.data.frame(dtn)$distance)
    d.gr$dist.feat[S4Vectors::queryHits(dtn)] = as.data.frame(dtn)$distance
  }
  cnv.df = data.frame(as.data.frame(ol.l), sample=cnv.gr$sample)
  ol.prof = unique(cnv.df[,1:length(feat.grl), drop=FALSE])
  null.gr = parallel::mclapply(1:nrow(ol.prof), function(ii){
    gr.ii = cnv.gr[which(apply(t(as.matrix(cnv.df[,1:length(feat.grl), drop=FALSE]))==unlist(ol.prof[ii,]),2,all))]
    widths = width(gr.ii)
    ## w.class = cut(widths, breaks=c(0,unique(quantile(widths, probs=ppoints(nb.class))), Inf))
    ## w.class = factor(w.class, levels=unique(w.class))
    w.uniq = unique(widths)
    if(length(w.uniq) < nb.class){
      w.class = widths
    } else {
      w.class = stats::cutree(stats::hclust(stats::dist(w.uniq), method="ward.D"), nb.class)
      names(w.class) = as.character(w.uniq)
      w.class = w.class[as.character(widths)]
    }
    tt = tapply(1:length(gr.ii), w.class, function(iii){
      w.i = widths[iii]
      demi.w.min = min(w.i, na.rm=TRUE) / 2
      demi.w.max = max(w.i, na.rm=TRUE) / 2
      good.d = colSums((t(as.matrix(d.df))-ifelse(unlist(ol.prof[ii,]),demi.w.min,demi.w.max))*ifelse(unlist(ol.prof[ii,]),1,-1)<0)
      nb.feat.good = ncol(ol.prof)
      while(all(good.d!=nb.feat.good)){
        nb.feat.good = nb.feat.good - 1
      }
      good.d = which(good.d==nb.feat.good)
      if(!is.null(dist.gr) & length(good.d) > length(w.i)){
        good.d = utils::head(good.d[sapply(gr.ii$dist.feat[iii], function(d)which.min(abs(d-d.gr$dist.feat[good.d])))], length(w.i))
      } else {
        good.d = sample(good.d, length(w.i), length(good.d) < length(w.i))
      }
      gr = GenomicRanges::resize(d.gr[good.d], w.i, fix="center")
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
  null.gr$dist.feat = NULL
  if(redo.duplicates && any((dup = duplicated(as.data.frame(null.gr))))){
    redo.gr = draw.controls(null.gr[which(dup)], feat.grl, nb.class=nb.class, nb.cores=nb.cores, redo.duplicates=FALSE, min.nb.gr=min.nb.gr, chr.prefix = chr.prefix, dist.gr=dist.gr)
    null.gr = c(null.gr[which(!dup)], redo.gr)
  }
  null.gr
    
}
