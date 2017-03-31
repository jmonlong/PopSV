##' Cluster samples using their abnormal regions. The similarity between two samples is computed as the amount of base pairs called in both samples divided by the average of the total base pairs called in each sample. The distance is 1 - this similarity.
##' @title Cluster samples using their abnormal regions
##' @param cnv.df a data.frame with the abnormal regions for each sample. Columns 'chr', 'start', 'end' and 'sample' are required.
##' @param cl.method the linkage criterion for hierarchical clustering. Default is 'complete'.
##' @param nb.cores the number of processors to use. Default is 1.
##' @return a list
##' \item{d}{a distance matrix with the distance between each pair of sample.}
##' \item{hc}{a 'hclust' object with the clustered samples.}
##' \item{cnv.geom}{if the geometric approximation was used, the geometric representation of each samples as a matrix.}
##' @author Jean Monlong
##' @export
cluster.regions <- function(cnv.df,cl.method="complete", nb.cores=3){
    cnv.gr = with(cnv.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), sample=sample))
    sample.names = unique(cnv.df$sample)
    dj.gr = GenomicRanges::disjoin(cnv.gr)
    cnv.geom = parallel::mclapply(sample.names, function(samp){
                                      dj.ol = IRanges::overlapsAny(dj.gr, cnv.gr[which(cnv.gr$sample==samp)])
                                      ifelse(dj.ol, GenomicRanges::width(dj.gr), 0)
                                  }, mc.cores=nb.cores)
    cnv.geom = matrix(unlist(cnv.geom), length(cnv.geom[[1]]))
    colnames(cnv.geom) = sample.names
    cnv.geom = cnv.geom / 1e3
    cnv.geom.sqrt = sqrt(cnv.geom)
    d = t(cnv.geom.sqrt) %*% cnv.geom.sqrt
    w.mean = apply(cnv.geom, 2, sum)
    sc = utils::combn(sample.names, 2)
    w.mean.mat = matrix(NA, length(w.mean), length(w.mean))
    w.mean.mat[lower.tri(w.mean.mat)] = w.mean[sc[1,]]+w.mean[sc[2,]]
    w.mean.mat[upper.tri(w.mean.mat)] = t(w.mean.mat)[upper.tri(w.mean.mat)]
    w.mean.mat = w.mean.mat / 2
    d.mat = 1 - d / w.mean.mat
    diag(d.mat) = 0
    hc = stats::hclust(stats::as.dist(d.mat),method=cl.method)
    return(list(d=d.mat,hc=hc, cnv.geom=cnv.geom))
}
