##' Compute P-values and Q-values from the Z-score distribution. 
##' @title P-values estimation
##' @param z a vector with the Z-scores
##' @param quant.int a vector with the quantiles to test for the estimation of the
##' null normal variance.
##' @return a list with
##' \item{pval}{the vector of P-values}
##' \item{qval}{the vector of Q-values / FDR estimates}
##' \item{quant.int}{the quantile used for the estimation of the null distribution
##' variance.}
##' \item{sigma.est}{the estimated null distribution variance}
##' @author Jean Monlong
##' @keywords internal
fdrtool.quantile <- function(z,quant.int = seq(.2,1,.02)){
    res = list(pval = rep(NA,length(z)),qval = rep(NA,length(z)), quant.int=NA, sigma.est=NA)
    z[which(is.infinite(z))] = NA ## Remove infinite values
    non.na.i = which(!is.na(z))
    z.non.na = z[non.na.i]
    sigma.est = sapply(quant.int, function(qi)fdrtool::censored.fit(z.non.na,quantile(abs(z.non.na),probs=qi,na.rm=TRUE))[5])
    res$quant.int = quant.int[which.min(abs(diff(sigma.est)))]
    res$sigma.est = sigma.est[which.min(abs(diff(sigma.est)))]
    pv = 2*pnorm(-abs(z.non.na),0,res$sigma.est)
    if(any(pv==0))
        pv[pv==0] = .Machine$double.xmin
    ft = fdrtool::fdrtool(pv,statistic="pvalue",plot=FALSE,verbose=FALSE)
    res$pval[non.na.i] = pv
    if(any(ft$qval<.05,na.rm=TRUE) | !any(pv<1e-10)){
        res$qval[non.na.i] = ft$qval
    } else {
        res$qval[non.na.i] = p.adjust(pv,method="fdr")
    }
    return(res)
}
