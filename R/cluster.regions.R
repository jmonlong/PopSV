##' Cluster samples with their abnormal regions. The similarity between two samples is computed as the amount of base pairs called in both samples divided by the average of the total base pairs called in each sample. The distance is 1 - this similarity.
##' @title Cluster samples with abnormal regions
##' @param cnv.df a data.frame with the abnormal regions for each sample. Columns 'chr', 'start', 'end' and 'sample' are required.
##' @param cl.method the linkage criterion for hierarchical clustering. Default is 'complete'.
##' @param nb.cores the number of processors to use. Default is 1.
##' @param geom.approx should the geometric approximation be used (faster). Default is FALSE.
##' @return a list
##' \item{d}{a distance matrix with the distance between each pair of sample.}
##' \item{hc}{a 'hclust' object with the clustered samples.}
##' \item{cnv.geom}{if the geometric approximation was used, the geometric representation of each samples as a matrix.}
##' @author Jean Monlong
##' @export
cluster.regions <- function(cnv.df,cl.method="complete", nb.cores=3, geom.approx=FALSE){
  cnv.gr = with(cnv.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), sample=sample))
  sample.names = unique(cnv.df$sample)
  if(geom.approx){
    dj.gr = GenomicRanges::disjoin(cnv.gr)
    cnv.geom = parallel::mclapply(sample.names, function(samp){
      dj.ol = GenomicRanges::overlapsAny(dj.gr, cnv.gr[which(cnv.gr$sample==samp)])
      ifelse(dj.ol, GenomicRanges::width(dj.gr), 0)
    }, mc.cores=nb.cores)
    cnv.geom = matrix(unlist(cnv.geom), length(cnv.geom[[1]]))
    colnames(cnv.geom) = sample.names
    cnv.geom = apply(cnv.geom, 2, function(x)x/sum(x))
    d.mat = as.matrix(dist(t(cnv.geom), method="man"))
  } else {
    dist.cnv <- function(samp1,samp2){
      cnv1.gr = cnv.gr[which(cnv.gr$sample==samp1)]
      cnv2.gr = cnv.gr[which(cnv.gr$sample==samp2)]
      ol = GenomicRanges::findOverlaps(cnv1.gr, cnv2.gr)
      w.ol = sum(GenomicRanges::width(GenomicRanges::pintersect(cnv1.gr[GenomicRanges::queryHits(ol)], cnv2.gr[GenomicRanges::subjectHits(ol)]))/1e3)
      1-(2 * w.ol / sum(GenomicRanges::width(c(cnv1.gr, cnv2.gr))/1e3))
    }
    sc = combn(sample.names, 2)
    d.l = parallel::mclapply(1:ncol(sc), function(ii) dist.cnv(sc[1,ii],sc[2,ii]), mc.cores=nb.cores)
    d.mat = matrix(NA,length(sample.names),length(sample.names))
    rownames(d.mat) = colnames(d.mat) = sample.names
    for(ii in 1:ncol(sc))
      d.mat[sc[1,ii],sc[2,ii]] = d.mat[sc[2,ii],sc[1,ii]] = d.l[[ii]]
    cnv.geom = NULL
  }
  hc = hclust(as.dist(d.mat),method=cl.method)
  return(list(d=d.mat,hc=hc, cnv.geom=cnv.geom))
}
