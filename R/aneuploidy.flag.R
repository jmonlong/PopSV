##' Chromosomes with clear aneuploidy are detected using a simple threshold-based approach. It is used to flag chromosome with complete aneuploidy that could handicap later analysis. However it might not detect partial chromosomal aberration or aneuploidy in samples with noisy read coverage.
##' @title Flag chromosomal aneuploidy
##' @param files.df a data.frame with information about samples such as the path to the count files.
##' @param col.file the name of the column in 'files.df' with the path information to use. Default is 'bc.gc.gz', i.e. GC corrected bin counts.
##' @param verbose Shoulf progress be outputted ? Default is TRUE.
##' @param ref.samples a vector with the names of the samples to be used as reference (e.g. normals). If NULL (default) all samples are used.
##' @param centromere.pos a vector with the position of the centromeres for each chromosome. Each position should be named with the chromosome name. If NULL (default), full chromosome aneuploidy is tested instead of arm-level deletion/duplication.
##' @param nb.cores the number of processing cores to use, Default is 1.
##' @return a data.frame with columns:
##' \item{sample}{the sample}
##' \item{chr}{the chromosome}
##' \item{pv.loss, pv.grain}{the p-value for a loss/gain}
##' \item{loss/gain}{is the loss/gain significant}
##' \item{bc.norm}{the normalized coverage}
##' \item{exp.bc}{the expected normalized coverage}
##' \item{bc.diff}{the difference between observed and expected coverage}
##' @author Jean Monlong
##' @export
aneuploidy.flag <- function(files.df, col.file = "bc.gc.gz", verbose=TRUE, ref.samples=NULL, centromere.pos=NULL, nb.cores=1) {
  if(is.null(ref.samples)){
      ref.samples = as.character(files.df$sample)
  }

  if(!is.null(centromere.pos) && is.null(names(centromere.pos))){
      stop("Positions in 'centromere.pos' should be named with the chromosome names")
  }

  if(all(colnames(files.df) != col.file)){
      stop("'col.file=", col.file, "' is not present in 'files.df'")
  }

  ## Compute median coverage per chromosome in each sample
  med.df = parallel::mclapply(1:nrow(files.df),function(file.ii){
                      if(verbose) message("Importing bin counts from sample ", file.ii, "...")
                      bc.df = utils::read.table(files.df[file.ii, col.file], as.is=TRUE, header=TRUE)
                      if(!is.null(centromere.pos)){
                          bc.df$arm = ifelse(bc.df$start<centromere.pos[as.character(bc.df$chr)], "p","q")
                          bc.df$chr = paste0(bc.df$chr, bc.df$arm)
                      }
                      data.frame(sample=files.df$sample[file.ii],
                                 stats::aggregate(bc~chr,bc.df, stats::median, na.rm=TRUE))
                  }, mc.cores=nb.cores)
  med.df = do.call(rbind, med.df)

  ## Adjust the median median counts
  if(verbose) message("Adjusting coverage...")
  med.med.df = stats::aggregate(bc~sample, med.df, stats::median, na.rm=TRUE)
  med.med.df$norm.fact = 1/med.med.df$bc
  med.med.df$bc = NULL
  med.df = merge(med.df, med.med.df)
  med.df$bc.norm = med.df$bc * med.df$norm.fact
  med.df$norm.fact = med.df$bc = NULL

  ## For each chr, fit null distribution on reference samples and estimate aneuploidy
  fit2norm.sd <- function(z, p0=c(p=.5,s1=1,s2=1)){ ## Fit 2 gaussian centered in 0
      mix.obj<-function(p,x){
          e<-p[1]*stats::dnorm(x/p[2])/p[2] + (1-p[1])*stats::dnorm(x/p[3])/p[3]
          if (any(e<=0) | p[1]<0 | p[1]>1) Inf
          else -sum(log(e))
      }
      lmix2a<-stats::deriv(~ -log(p*stats::dnorm(x/s1)/s1 + (1-p)*stats::dnorm(x/s2)/s2), c("p","s1","s2"), function(x,p,s1,s2) NULL)
      mix.gr<-function(pa,x){
          p<-pa[1]
          s1<-pa[2]
          s2<-pa[3]
          colSums(attr(lmix2a(x,p,s1,s2),"gradient"))}
      results=stats::optim(p0,mix.obj,mix.gr,x=z)
      data.frame(p=results$par[1], s1=results$par[2], s2=results$par[3])
  }
  aneuChr <- function(df){ ## analyze one chromosome
      df$exp.bc = stats::median(df$bc.norm[which(df$sample %in% ref.samples)], na.rm=TRUE)
      df$bc.diff = df$bc.norm - df$exp.bc
      fit.res = fit2norm.sd(df$bc.diff[which(df$sample %in% ref.samples)])
      df$pv.loss = stats::pnorm(df$bc.diff, 0, fit.res$s1)
      df$pv.gain = 1-stats::pnorm(df$bc.diff, 0, fit.res$s1)
      df$loss = stats::p.adjust(df$pv.loss)<.01
      df$gain = stats::p.adjust(df$pv.gain)<.01
      df
  }
  if(verbose) message("Testing for aneuploidy...")
  med.fit = parallel::mclapply(unique(med.df$chr), function(chr.in)aneuChr(med.df[which(med.df$chr == chr.in),]), mc.cores=nb.cores)
  med.fit = do.call(rbind, med.fit)

  return(med.fit)
}
