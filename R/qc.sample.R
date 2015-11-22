##' Compute QC metrics for each sample informing how much it fits the reference samples.
##' @title QC samples
##' @param bins.df a data.frame with information about the bins. Columns 'chr', 'start', 'end' are required.
##' @param files.df a data.frame with information about the files, usually created by 'init.filenames'.
##' @param cnv.df a data.frame with PopSV calls, created by 'call.abnormal.cov'.
##' @param ref.samples a vector with the names of the samples used as reference. If NULL (default), all samples in 'files.df' are used.
##' @param n.subset the number of bins to use to compute the metrics. Default is 1000.
##' @param nb.cores number of cores to use. If higher than 1, \code{parallel}
##' package is used to parallelized the counting.
##' @return a data.frame with qc metrics for each sample of the input 'files.df'.
##' @author Jean Monlong
##' @export
qc.sample <- function(bins.df, files.df=NULL, cnv.df=NULL, ref.samples=NULL, n.subset=1000, nb.cores=1){

  if(!all(c("chr","start","end") %in% colnames(bins.df))){
    stop("Missing columns in 'bins.df'. Columns 'chr', 'start', 'end' are required.")
  }
  if(is.null(files.df) & is.null(cnv.df)){
    stop("At least one of 'files.df' or 'cnv.df' parameters must be non-NULL.")
  }
  if(is.null(ref.samples) & !is.null(files.df)){
    ref.samples = files.df$sample
  }


  res = NULL
  if(!is.null(files.df)){
    bc.df = quick.count(files.df, bins.df, col.files="bc.gc.gz", nb.rand.bins=n.subset, nb.cores=nb.cores)
    samples = files.df$sample
    bc.df = med.norm(bc.df, nb.cores=nb.cores, norm.stats.comp=FALSE)$bc.norm
    d.geom = as.matrix(dist(t(bc.df[,samples])))
    colnames(d.geom) = samples
    ## d.o = rowMeans(d.geom[,ref.samples])/mean(d.geom[,ref.samples])
    d.o = apply(d.geom[,ref.samples], 1, function(d.f) mean(d.f[d.f<quantile(d.f, probs=.1, na.rm=TRUE)], na.rm=TRUE)) / mean(d.geom[d.geom<quantile(d.geom, probs=.1, na.rm=TRUE)], na.rm=TRUE)
    res = data.frame(sample=samples, dist.bc.ref=d.o)

    if(any(colnames(bins.df)=="m")){
      lowcov.bins.df = bins.df[order(bins.df$m),]
      lowcov.bins.df = lowcov.bins.df[head(which(!is.na(lowcov.bins.df$sd)),n.subset),]
      bc.df = quick.count(files.df, lowcov.bins.df, col.files="bc.gc.gz", nb.rand.bins=n.subset, nb.cores=nb.cores)
      d.geom = as.matrix(dist(scale(t(scale(bc.df[,samples], scale=FALSE)))))
      colnames(d.geom) = samples
      d.o = rowMeans(d.geom[,ref.samples])/mean(d.geom[,ref.samples])
      res$dist.bc.ref.lowcov = d.o
    }
  }

  if(!is.null(cnv.df)){
    nb.bin.cons = fc = NULL ## Uglily silence R checks
    res.cnv = dplyr::summarize(dplyr::group_by(cnv.df, sample), single.bin.prop=mean(nb.bin.cons==1, na.rm=TRUE), cn2.prop = mean(abs(fc[which(nb.bin.cons>2)])<.25, na.rm=TRUE))
    if(!is.null(res)){
      res = merge(res, res.cnv, all=TRUE)
      rownames(res) = NULL
    } else {
      res = res.cnv
    }
  }

  return(res)
}
