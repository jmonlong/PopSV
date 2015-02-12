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
z.comp <- function(files.df, samples, z.poisson=FALSE, col="bc.gc.tnk"){
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
  bc.df = createEmptyDF(c("character",rep("integer",2), rep("numeric",length(samples))), nrow(bc.1))
  z.df = createEmptyDF(c("character",rep("integer",2), rep("numeric",length(samples))), nrow(bc.1))
  fc.df = createEmptyDF(c("character",rep("integer",2), rep("numeric",length(samples))), nrow(bc.1))
  msd.df = createEmptyDF(c("character",rep("integer",2), rep("numeric",3)), nrow(bc.1))
  colnames(bc.df) = colnames(z.df) = colnames(fc.df) = c("chr","start","end", samples)
  colnames(msd.df) = c("chr","start","end", "m","sd","nb.remove")
  bc.df$chr = z.df$chr = fc.df$chr = bc.df$chr
  bc.df$start = z.df$start = fc.df$start = bc.df$start
  bc.df$end = z.df$end = fc.df$end = bc.df$end
  bc.df[,samples[1]] = bc.1[,4]
  if(nb.cores>1){
    bc.l = parallel::mclapply(subset(files.df, sample%in%samples[-1])[,col.bc], function(fi){
      read.table(fi,header=TRUE, as.is=TRUE)[,4]
    },mc.cores=nb.cores)
  } else {
    bc.l = lapply(subset(files.df, sample%in%samples[-1])[,col.bc], function(fi){
      read.table(fi,header=TRUE, as.is=TRUE)[,4]
    })
  }
  for(samp.i in 1:(length(samples)-1)){
    bc.df[,as.character(samples[1+samp.i])] = bc.l[[samp.i]]
  }
  for(ii in 1:nrow(bc.df)){
    bc.t = bc.df[ii,-(1:3)]
    msd = mean.sd.outlierR(bc.t,1e-6)
    msd.df[ii, -(1:3)] = with(msd, c(m, sd, nb.remove))
    z.df[ii,-(1:3)] = z.comp.f(bc.t,msd$m,msd$sd)
    fc.df[ii,-(1:3)] = bc.t/msd$m
  }
  return(list(z=z.df, fc=fc.df, msd=msd.df, z.poisson=z.poisson))
}
