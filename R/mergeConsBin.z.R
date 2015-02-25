##' Merge consecutive abnormal bins into one larger event. Two consecutive bins
##' will be merged if the maximum (minimum) of their two Z-scores if significantly
##' lower (higher) than what is observed from a null distribution. The null distribution
##' is simulated by taking the maximum (minimum) values of two Normal values. The variance
##' of the simulated Normal is given as a parameter of the function.
##' @title Merge abnormal consecutive bins
##' @param res.df a data.frame with the Z-scores. Columns 'chr', 'start', 'end' and 'z'
##' are required.
##' @param fdr.th the False Discovery Rate threshold.
##' @param col.mean a list of potential columns for which the mean value should be kept
##' when merging bins.
##' @param sd.null the estimated standard deviation of the Z-score null distribution.
##' Usually, computed during P-value estimation by 'fdrtool.quantile'.
##' @param nb.sim the number of simulated Z-scores for the P-value computation.
##' @return a data.frame similar to the input 'res.df' but with an extra 'nb.bin.cons'
##' column (the number of bin merged for each event).
##' @author Jean Monlong
##' @keywords internal
mergeConsBin.z <- function(res.df,fdr.th=.05,col.mean=c("z","pv","qv","fc"),sd.null=1,nb.sim=1e6){
    ## Compute P-value from an empirical null distribution
    compute.pv <- function(z,z.null, alt.greater=TRUE){
        if(length(z)>1){
            z.f = factor(z)
            zn.c = cut(z.null,c(-Inf,levels(z.f),Inf),right=FALSE)
            zn.sum = table(zn.c)
            zn.cs = cumsum(zn.sum)[-length(zn.sum)]
            names(zn.cs) = levels(z.f)
            zn.cs.f = zn.cs[z.f] ## Pb ? as.character, luckily not
            names(zn.cs.f) = names(z)
            if(alt.greater){
                return(1 - ( zn.cs.f / (length(z.null)+1) ))
            } else {
                return((zn.cs.f+1) / (length(z.null)+1))
            }
        } else {
            if(alt.greater){
                pv = (sum(z <= z.null) + 1) / (length(z.null) + 1)
            } else {
                pv = (sum(z >= z.null) + 1) / (length(z.null) + 1)
            }
            names(pv) = names(z)
            return(pv)
        }
    }
    ## Simulate the empirical null distribution
    z.null = apply(rbind(rnorm(nb.sim,0,sd.null),rnorm(nb.sim,0,sd.null)),2,sort)
    ## Compute P-values
    res.df = with(res.df, dplyr::arrange(res.df, chr, start))
    pvLink.f <- function(df){
        z.link = apply(rbind(df$z[-1],df$z[-nrow(df)]),2,sort)
        data.frame(pv.dup = compute.pv(z.link[1,],z.null=z.null[1,]),
                   pv.del = compute.pv(z.link[2,],z.null=z.null[2,],alt.greater=FALSE))
    }
    chr = . = link = NULL ## Uglily appease R checks
    link.df = dplyr::do(dplyr::group_by(res.df,chr),pvLink.f(.))
    ## Multiple test correction
    link.df$qv.dup = fdrtool::fdrtool(link.df$pv.dup, statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
    link.df$qv.del = fdrtool::fdrtool(link.df$pv.del, statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
    ## Annotate and merge bins
    link.annotate.f <- function(df){
      cat(df$chr[1],"\n")
      link.df = link.df[which(link.df$chr==df$chr[1]),]
        link.v = rep("none", nrow(df))
        link.v[c(FALSE,link.df$qv.dup<=fdr.th) | c(link.df$qv.dup<=fdr.th,FALSE) | (df$qv<=fdr.th & df$z>0)] = "dup"
        link.v[c(FALSE,link.df$qv.del<=fdr.th) | c(link.df$qv.del<=fdr.th,FALSE) | (df$qv<=fdr.th & df$z<0)] = "del"
        rl = rle(link.v)
        rlv = rl$values
        rlv[rlv!="none"] = paste(link.df$chr[1],rlv[rlv!="none"],1:sum(rlv!="none"))
        df$link = rep(rlv, rl$lengths)
        return(df)
    }
    res.df = dplyr::do(dplyr::group_by(res.df,chr),link.annotate.f(.))
    res.df = res.df[which(res.df$link!="none"),]
    if(nrow(res.df)>0){
        merge.bin.f <- function(df){
            res = data.frame(chr=df$chr[1],
                start=df$start[1],
                end=df$end[nrow(df)],
                nb.bin.cons=nrow(df))
            if(nrow(df)>2) df = df[c(2,nrow(df)-1),]
            cbind(res,t(apply(df[,intersect(colnames(df),col.mean),drop=FALSE],2,mean)))
        }
        res.df = dplyr::do(dplyr::group_by(res.df,link),merge.bin.f(.))
    } else {
        res.df = NULL
    }
    res.df$link = NULL
    return(as.data.frame(res.df))
}
