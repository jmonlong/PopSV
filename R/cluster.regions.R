##' Cluster samples with their abnormal regions.
##' @title Cluster samples with abnormal regions
##' @param cnv.df a data.frame with the abnormal regions for each sample. Columns 'chr', 'start', 'end' and 'sample' are required.
##' @param cl.method the linkage criterion for hierarchical clustering. Default is 'complete'.
##' @param nb.cores the number of processors to use. Default is 1. 
##' @return a list
##' \item{d}{a distance matrix with the distance between each pair of sample.}
##' \item{hc}{a 'hclust' object with the clustered samples.}
##' @author Jean Monlong
##' @export
cluster.regions <- function(cnv.df,cl.method="complete", nb.cores=3){
  cnv.gr = with(cnv.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), sample=sample))
  dist.cnv <- function(samp1,samp2){
    cnv1.gr = cnv.gr[which(cnv.gr$sample==samp1)]
    cnv2.gr = cnv.gr[which(cnv.gr$sample==samp2)]
    ol = GenomicRanges::findOverlaps(cnv1.gr, cnv2.gr)
    w.ol = sum(GenomicRanges::width(GenomicRanges::pintersect(cnv1.gr[GenomicRanges::queryHits(ol)], cnv2.gr[GenomicRanges::subjectHits(ol)]))/1e3)
    1-(2 * w.ol / sum(GenomicRanges::width(c(cnv1.gr, cnv2.gr))/1e3))
  }
  sample.names = unique(cnv.df$sample)
  sc = combn(sample.names, 2)
  d.l = parallel::mclapply(1:ncol(sc), function(ii) dist.cnv(sc[1,ii],sc[2,ii]), mc.cores=nb.cores)
  d.mat = matrix(NA,length(sample.names),length(sample.names))
  rownames(d.mat) = colnames(d.mat) = sample.names
  for(ii in 1:ncol(sc))
    d.mat[sc[1,ii],sc[2,ii]] = d.mat[sc[2,ii],sc[1,ii]] = d.l[[ii]]
  hc = hclust(as.dist(d.mat),method=cl.method)
  return(list(d=d.mat,hc=hc))
}
