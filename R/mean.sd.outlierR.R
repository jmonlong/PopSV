##' Computes the mean and standard deviation after removing outlier groups.
##' \code{outliers} package and Grubbs approach are used.
##' @title Mean and standard diviation after outlier removal
##' @param x a vector of numeric values.
##' @param pv.max.ol the maximum P-value for the Grubbs test. Any outlier
##' with a P-value lower will be considered as outlier.
##' @return a vector similar to 'x' but with NAs for outliers.
##' @author Jean Monlong
##' @keywords internal
mean.sd.outlierR <- function(x,pv.max.ol=1e-6){
    rm.outOl <- function(x,pv.max=1e-6){
        grubbs.sw <- function(x,type=10,opposite=FALSE,two.sided=TRUE,max.pv=.01,max.step=10){
            grubbs.t <- function(x,type=10,opposite=FALSE,two.sided=TRUE){
                gt = outliers::grubbs.test(x,type,opposite,two.sided)
                gt$outlier = outliers::outlier(x)
                gt$outlier.i = which(outliers::outlier(x,logical=TRUE))[1]
                return(gt)
            }
            pv = outliers = NULL
            step = 1
            continue = TRUE
            while(continue & step <= max.step){
                gt = grubbs.t(as.numeric(x),type,opposite,two.sided)
                if(gt$p.value>max.pv){
                    continue=FALSE
                } else {
                    pv = c(pv, gt$p.value)
                    outliers = c(outliers, gt$outlier.i)
                    x[outliers[step]] = NA
                }
                step = step + 1
            }
            return(list(x=x,pv=pv,outliers=outliers,max.pv=max.pv))
        }
        xs = sort(x)
        dd = diff(xs)
        gt = grubbs.sw(dd,max.pv=pv.max)
        if(!is.null(gt$outliers)){
            cut.pos = c(0,sort(gt$outliers),length(xs))
            gp.sizes = diff(cut.pos)
            cm = which.max(gp.sizes)
            xss = rep(NA,length(xs))
            xss[(cut.pos[cm]+1):cut.pos[cm+1]] = xs[(cut.pos[cm]+1):cut.pos[cm+1]]
        } else {
            xss = xs
        }
        xss
    }
    sd.mad <- function(x){
        if(all(x==0,na.rm=TRUE)) x = rpois(length(x),.5)
        sd.res = mad(x, na.rm=TRUE)
        if(sd.res==0) return(sd.mad(rep(0,length(x))))
        else return(sd.res)
    }
    trim.mean <- function(x){
        x.ii = order(abs(x-median(x,na.rm=TRUE)))[1:10]
        return(mean(x[x.ii], na.rm=TRUE))
    }
    x.rm = rm.outOl(x,pv.max=pv.max.ol)
    return(list(m=trim.mean(x.rm),sd=sd.mad(x.rm),nb.remove=length(x)-sum(!is.na(x.rm))))
}
