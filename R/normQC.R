##' Quality control statistics of bin counts. First, the counts are tested to follow a
##' normal (or Poisson) distribution across samples, in each bin. Then, the randomness
##' of sample ranks are tested. Finally, Z-scores are computed for a subset of the bins
##' and their normality is tested. The second and third test are the most important.
##' Indeed, consistent rankings supports sample-specific technical bias, hence reduced
##' power to detect "true" abnormal read counts. Non-normal Z-scores will lead to
##' inappropriate fit for the null distribution.
##'
##' Shapiro test is used to test normality of the bin counts across samples. The proportion
##' of bins with non-normal distribution is derived from the Pi0 estimate estimated by
##' package \code{qvalue}. Pi0 is the proportion of pvalues following the null distribution.
##'
##' Goodness of fir from package \code{vcd} is used to test if bin counts follow a Poisson
##' distribution. Again Pi0 estimate from \code{qvalue} package is used to compute the proportion
##' of bin that don't follow Poisson distribution.
##'
##' The randomness of the sample ranks are tested using a Chi-Square test. To have an
##' idea of the proportion of the bins with abnormal ranking, bins are randomly divided
##' into groups and the test performed in each group. The proportion of bin with non-random
##' ranks is approximated by the proportion of groups which failed the Chi-Square test.
##'
##' Z-scores normality is computed by comparing their density distribution and a fitted
##' normal distribution. The estimate represents the proportion of the area under the curve
##' that is unique to the Z-score curve.
##' @title Normalized bin count QC metrics
##' @return a list with
##' \item{prop.non.normal.bin}{proportion of bins with non-normal distribution across samples.}
##' \item{prop.non.poisson.bin}{proportion of bins with non-Poisson distribution across samples.}
##' \item{prop.nonRand.rank}{proportion of bins with non-random ranks.}
##' \item{prop.non.norm.z.mean}{average (across samples) proportion of bins with non-random Z-scores.}
##' \item{prop.non.norm.z.max}{maximum (i.e for worst sample) proportion of bins with non-random Z-scores.}
##' \item{prop.non.norm.z}{proportion of bins with non-random Z-scores, for each sample.}
##' \item{n.subset}{number of bins used for the analysis.}
##' @param bc.mat matrix with bin counts (bins x samples). 
##' @param n.subset number of bins to use for the analysis. Default is 10 000. Bins are selected randomly.
##' @author Jean Monlong
##' @export
normQC <- function(bc.mat, n.subset=1e4){
    ## library(MASS)
    ## library(vcd)
    ## library(qvalue)
    
  if(!is.null(n.subset) & n.subset<nrow(bc.mat)){
    bc.mat = bc.mat[sample(1:nrow(bc.mat), n.subset),]
  }
  
  ## Bin count normality
  pv.normal = apply(bc.mat, 1, function(bc.i){
    if(length(unique(bc.i))==1) return(1)
    return(shapiro.test(bc.i)$p.value)
  })
  pv.normal = as.numeric(na.omit(pv.normal))
  qv.normal = qvalue::qvalue(pv.normal)

  ## Bin count Poisson
  pv.poisson = apply(bc.mat, 1, function(bc.i){
    return(as.numeric(summary(vcd::goodfit(bc.i, type="poisson"))[1,3]))
  })
  pv.poisson = as.numeric(na.omit(pv.poisson))
  qv.poisson = qvalue::qvalue(pv.poisson)

  ## Ranks randomness
  win.size = 100
  rank.mat = apply(bc.mat, 1, function(bc.i){
    if(any(!is.na(bc.i))){
      rank(bc.i,ties.method="random")
    } else {
      rep(NA,length(bc.i))
    }
  })
  group.bins = rep(1:ceiling(nrow(bc.mat)/win.size),each=win.size)[1:nrow(bc.mat)]
  pv.rank = tapply(1:nrow(bc.mat), group.bins, function(bin.id){
    rk = rank.mat[,bin.id]
    if(length((noNA = which(!is.na(rk[1,]))))>1){
      rk = rk[,noNA] ## Remove NA
      rk.t = apply(rk,1,function(e)table(factor(e,levels=1:nrow(rk))))
      return(suppressWarnings(chisq.test(rk.t)$p.value))
    } else {
      return(NA)
    }
  })
  
  ## Z-score distribution
  n.dens = 1000
  z = apply(bc.mat,1,function(bc.i){
    if(any(is.na(bc.i))){
      return(rep(NA, length(bc.i)))
    }
    bc.i <- as.numeric(bc.i)
    (bc.i-mean(bc.i,na.rm=TRUE))/sd(bc.i,na.rm=TRUE)
  })
  rownames(z) = colnames(bc.mat)
  non.norm.z = apply(z, 1, function(z.samp){
    z.samp = as.numeric(na.omit(z.samp))
    z.norm.est = MASS::fitdistr(z.samp,"normal")$estimate
    f = density(z.samp,n=n.dens,from=min(z.samp),to=max(z.samp))
    fn = dnorm(seq(min(z.samp),max(z.samp),length.out=n.dens),z.norm.est[1],z.norm.est[2])
    dis.prop = sum(abs(f$y-fn)) / (2 * sum(fn))
    return(dis.prop)
  })
  
  ## 
  return(list(prop.non.normal.bin = 1-qv.normal$pi0,
              prop.non.poisson.bin = 1-qv.poisson$pi0,
              prop.nonRand.rank = mean(pv.rank<=.05),
              prop.non.norm.z.mean = mean(non.norm.z),
              prop.non.norm.z.max = max(non.norm.z),
              prop.non.norm.z = non.norm.z,
              n.subset = n.subset))
}
