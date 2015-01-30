##' Compute P-values and Q-values from the Z-score distribution. 
##' @title P-values estimation from mixture of 2 centered normal
##' @param z a vector with the Z-scores
##' @param plot should some graphs be displayed. Default if FALSE.
##' @param min.prop.null the minimum proportion of the genome expected to be normal.Recommended: 0.8 for "normal" samples, 0.4 (default) for tumor samples.
##' @return a list with
##' \item{pval}{the vector of P-values}
##' \item{qval}{the vector of Q-values / FDR estimates}
##' \item{sigma.est.dup}{the estimated null distribution variance for positive Z-scores}
##' \item{sigma.est.del}{the estimated null distribution variance for negative Z-scores}
##' @author Jean Monlong
##' @keywords internal
fdrtool.quantile.2N <- function(z, plot=TRUE, min.prop.null=.95){

  localMax <- function(x,min.max.prop=.1){
    d = density(x,na.rm=TRUE)
    im = 1+which(diff(sign(diff(d$y)))==-2)
    my = max(d$y)
    max.id = im[which(d$y[im] >= min.max.prop * my)]
    max.id.o = max.id[order(d$y[max.id],decreasing=TRUE)]
    return(list(lM=d$x[max.id.o], h=d$y[max.id.o]/my))
  }
  fit2norm.sd.cens <- function(z, p0=c(p=.75,s1=1,s2=2), z0){
    z = z[abs(z)<z0]
    mix.obj<-function(p,x){
      e<-p[1]*dnorm(x/p[2])/((pnorm(z0,0,p[2])-pnorm(-z0,0,p[2]))*p[2]) + (1-p[1])*dnorm(x/p[3])/((pnorm(z0,0,p[3])-pnorm(-z0,0,p[3]))*p[3])
      if (any(e<=0, na.rm=TRUE) | p[1]<0 | p[1]>1) Inf
      else -sum(log(e))
    }
    lmix2a<-deriv(~ -log(p*dnorm(x/s1)/((pnorm(z0,0,s1)-pnorm(-z0,0,s1))*s1) + (1-p)*dnorm(x/s2)/((pnorm(z0,0,s2)-pnorm(-z0,0,s2))*s2)), c("p","s1","s2"), function(x,p,s1,s2) NULL)
    mix.gr<-function(pa,x){
      p<-pa[1]
      s1<-pa[2]
      s2<-pa[3]
      colSums(attr(lmix2a(x,p,s1,s2),"gradient"))}
    results=optim(p0,mix.obj,mix.gr,x=z)
    if(results$par[1]<.5){
      results$par[1] = 1-results$par[1]
      results$par[2:3] = results$par[3:2]
    }
    results
  }
  sim2norm.sd <- function(pars, nb.sims=1e6){
    if(pars["p"]<0 | pars["p"]>1)
      return(rnorm(nb.sims, 0, pars["s1"]))
    c(rnorm(pars["p"]*nb.sims, 0, pars["s1"]),
      rnorm((1-pars["p"])*nb.sims, 0, pars["s2"]))
  }
  find.par <- function(z, z0.int=seq(quantile(abs(z), probs=min.prop.null),quantile(abs(z), probs=.99),.2)){
    p.df = plyr::ldply(z0.int, function(z0){
      p = fit2norm.sd.cens(z, z0=z0)
      zsim = sim2norm.sd(p$par)
      dz = density(z[abs(z)<z0], from=-z0, to=z0, n=512)
      dzsim = density(zsim[abs(zsim)<z0], from=-z0, to=z0, n=512)
      data.frame(z0=z0, dens.diff = sum(abs(dz$y-dzsim$y))/sum(dzsim$y), p=p$par[1], s1=p$par[2], s2=p$par[3])
    })
    dd = localMax(p.df$dens.diff)
    list(par=unlist(p.df[which.min(abs(p.df$dens.diff-min(dd$lM))),]))
  }
  p2norm <- function(z, pars){
    pars["p"]*pnorm(z,0,pars["s1"]) + (1-pars["p"])*pnorm(z,0,pars["s2"])
  }
  
  res = list(pval = rep(NA,length(z)),qval = rep(NA,length(z)), sigma.est.dup=NA, sigma.est.del=NA)
  z[which(is.infinite(z))] = NA ## Remove infinite values
  non.na.i = which(!is.na(z) & z!=0)
  z.non.na = z[non.na.i]

  sup.ss = 5e4
  ## Duplication
  z.dup = z.non.na[z.non.na>0]
  z.dup = sample(c(-1,1),length(z.dup), replace=TRUE)*z.dup
  if(length(z.dup)>sup.ss){
    p = find.par(sample(z.dup, sup.ss))
  } else {
    p = find.par(z.dup)   
  }
  res$pval[non.na.i[z.non.na>0]] = 2*p2norm(-abs(z.dup), p$par)
  ## Deletion
  z.del = z.non.na[z.non.na<0]
  z.del = sample(c(-1,1),length(z.del), replace=TRUE)*z.del
  if(length(z.del)>sup.ss){
    p = find.par(sample(z.del, sup.ss))
  } else {
    p = find.par(z.del)
  }
  res$pval[non.na.i[z.non.na<0]] = 2*p2norm(-abs(z.del), p$par)

  if(any(res$pval==0, na.rm=TRUE))
    res$pval[which(res$pval==0)] = .Machine$double.xmin
  res$qval = p.adjust(res$pval,method="fdr")

  if(plot & any(!is.na(res$pval))){
    plot.df = data.frame(z=z, pv=res$pval, qv=res$qval)
    print(ggplot2::ggplot(subset(plot.df,abs(z)<quantile(abs(z), probs=.95)+1),ggplot2::aes(x=z)) +
          ggplot2::geom_histogram() + 
          ggplot2::xlab("Z-score") + 
          ggplot2::ylab("number of bins") + 
          ggplot2::theme_bw())
    print(ggplot2::ggplot(plot.df,ggplot2::aes(x=pv)) + ggplot2::geom_histogram() +
          ggplot2::xlab("P-value") + ggplot2::xlim(0,1) + 
          ggplot2::ylab("number of bins") + 
          ggplot2::theme_bw())
    print(ggplot2::ggplot(subset(plot.df,qv<.1),ggplot2::aes(x=cut(qv, breaks=c(-Inf,.001,.01,.5,.1)))) + ggplot2::geom_bar() +
          ggplot2::xlab("Q-value") + 
          ggplot2::ylab("number of bins") + 
          ggplot2::theme_bw())
  }
  
  return(res)
}
