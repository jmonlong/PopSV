##' @title Z-score computation
##' @return a list with
##' \item{z}{a data.frame with the Z-scores for each bin and sample (bin x sample).}
##' \item{fc}{a data.frame with the fold-change compared to the average bin count in
##' the reference samples for each bin and sample (bin x sample).}
##' \item{msd}{the mean, standard deviation and number of removed outlier samples in each bin.}
##' @author Jean Monlong
##' @export
##' @param files.df a data.frame with the file paths.
##' @param samples the samples to analyze.
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @param col the column in 'files.df' that define the bin count file path.
##' @param nb.cores the number of cores to use.
z.comp <- function(files.df, samples, z.poisson=FALSE, col="bc.gc.tnk", nb.cores=1){
  if(z.poisson){
    z.comp.f <- function(x, mean.c, sd.c){
      z.n = (x-mean.c)/sd.c
      z.p = qnorm(ppois(x, mean.c))
      n.ii = abs(z.n) < abs(z.p)
      z.p[n.ii] = z.n[n.ii]
      z.p
    }
  } else {
    z.comp.f <- function(x, mean.c, sd.c){(x-mean.c)/sd.c}
  }
  bc.1 = read.table(subset(files.df, sample==samples[1])[,col],header=TRUE, as.is=TRUE)
  if(nb.cores>1){
    bc.l = parallel::mclapply(subset(files.df, sample%in%samples)[,col], function(fi){
      read.table(fi,header=TRUE, as.is=TRUE)[,4]
    },mc.cores=nb.cores)
  } else {
    bc.l = lapply(subset(files.df, sample%in%samples)[,col], function(fi){
      read.table(fi,header=TRUE, as.is=TRUE)[,4]
    })
  }
  bc.l = matrix(unlist(bc.l), length(bc.l[[1]]))
  msd = apply(bc.l, 1, function(ee)unlist(mean.sd.outlierR(ee, pv.max.ol=1e-6)))
  z = apply(bc.l, 2, z.comp.f, mean.c=msd[1,], sd.c=msd[2,])
  fc = bc.l/msd[1,]
  colnames(z) = colnames(fc) = samples
  z = data.frame(bc.1[,1:3],z)
  fc = data.frame(bc.1[,1:3],fc)
  msd = data.frame(bc.1[,1:3],t(msd))
  return(list(z=z, fc=fc, msd=msd, z.poisson=z.poisson))
}
