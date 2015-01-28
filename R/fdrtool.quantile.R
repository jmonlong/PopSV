##' Compute P-values and Q-values from the Z-score distribution. 
##' @title P-values estimation
##' @param z a vector with the Z-scores
##' @param quant.int a vector with the quantiles to test for the estimation of the
##' null normal variance.
##' @param ref.dist.weight the weight (value between 0 and 1) based on the distance to the reference samples.
##' @return a list with
##' \item{pval}{the vector of P-values}
##' \item{qval}{the vector of Q-values / FDR estimates}
##' \item{sigma.est.dup}{the estimated null distribution variance for positive Z-scores}
##' \item{sigma.est.del}{the estimated null distribution variance for negative Z-scores}
##' @author Jean Monlong
##' @keywords internal
fdrtool.quantile <- function(z,quant.int = seq(.4,1,.02), ref.dist.weight=NULL){
    localMax <- function(x,min.max.prop=.1){
        d = density(x,na.rm=TRUE)
        im = 1+which(diff(sign(diff(d$y)))==-2)
        my = max(d$y)
        max.id = im[which(d$y[im] >= min.max.prop * my)]
        max.id.o = max.id[order(d$y[max.id],decreasing=TRUE)]
        return(list(lM=d$x[max.id.o], h=d$y[max.id.o]/my))
    }
    
    res = list(pval = rep(NA,length(z)),qval = rep(NA,length(z)), sigma.est.dup=NA, sigma.est.del=NA)
    z[which(is.infinite(z))] = NA ## Remove infinite values
    non.na.i = which(!is.na(z) & z!=0)
    z.non.na = z[non.na.i]

    ## Duplication
    z.dup = z.non.na[z.non.na>0]
    z.dup = sample(c(-1,1),length(z.dup), replace=TRUE)*z.dup
    sd.df = plyr::ldply(quant.int, function(qi){
      data.frame(quant=qi, sd.est = fdrtool::censored.fit(z.dup,quantile(abs(z.dup),probs=qi,na.rm=TRUE))[5])
    })
    sd.e = localMax(sd.df$sd.est)$lM[1]
    if(!is.null(ref.dist.weight)){
      sd.df = dplyr::arrange(sd.df, abs(sd.est-sd.e))
      sd.e = sd.e + ref.dist.weight*10*sd(sd.df$sd.est[1:20])
    }
    res$sigma.est.dup = sd.e
    ## Deletion
    z.del = z.non.na[z.non.na<0]
    z.del = sample(c(-1,1),length(z.del), replace=TRUE)*z.del
    sd.df = plyr::ldply(quant.int, function(qi){
      data.frame(quant=qi, sd.est = fdrtool::censored.fit(z.del,quantile(abs(z.del),probs=qi,na.rm=TRUE))[5])
    })
    sd.e = localMax(sd.df$sd.est)$lM[1]
    if(!is.null(ref.dist.weight)){
      sd.df = dplyr::arrange(sd.df, abs(sd.est-sd.e))
      sd.e = sd.e + ref.dist.weight*10*sd(sd.df$sd.est[1:20])
    }
    res$sigma.est.del = sd.e

    res$pval[non.na.i[z.non.na>0]] = 2*pnorm(-z.dup,0,res$sigma.est.dup)
    res$pval[non.na.i[z.non.na<0]] = 2*pnorm(z.del,0,res$sigma.est.del)
    if(any(res$pval==0, na.rm=TRUE))
        res$pval[which(res$pval==0)] = .Machine$double.xmin
    res$qval = p.adjust(res$pval,method="fdr")
    return(res)
}
