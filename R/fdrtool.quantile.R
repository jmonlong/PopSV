##' Compute P-values and Q-values from the Z-score distribution. 
##' @title P-values estimation
##' @param z a vector with the Z-scores
##' @param quant.int a vector with the quantiles to test for the estimation of the
##' null normal variance.
##' @param ref.dist.weight the weight (value between 0 and 1) based on the distance to the reference samples.
##' @param plot should some graphs be displayed. Default if FALSE.
##' @param z.th the Z-score threshold to use to estimate the null normal variance, if 'quant.int=NULL'. If one value the same threshold is used for both duplication and deletion; if a vector, the threshold for duplication / deletions. 
##' @return a list with
##' \item{pval}{the vector of P-values}
##' \item{qval}{the vector of Q-values / FDR estimates}
##' \item{sigma.est.dup}{the estimated null distribution variance for positive Z-scores}
##' \item{sigma.est.del}{the estimated null distribution variance for negative Z-scores}
##' @author Jean Monlong
##' @keywords internal
fdrtool.quantile <- function(z, quant.int = seq(0.4, 1, 0.02), ref.dist.weight = NULL, 
    plot = TRUE, z.th=NULL) {
    localMax <- function(x, min.max.prop = 0.1) {
        d = density(x, na.rm = TRUE)
        im = 1 + which(diff(sign(diff(d$y))) == -2)
        my = max(d$y)
        max.id = im[which(d$y[im] >= min.max.prop * my)]
        max.id.o = max.id[order(d$y[max.id], decreasing = TRUE)]
        return(list(lM = d$x[max.id.o], h = d$y[max.id.o]/my))
    }
    
    res = list(pval = rep(NA, length(z)), qval = rep(NA, length(z)), sigma.est.dup = NA, 
        sigma.est.del = NA)
    z[which(is.infinite(z))] = NA  ## Remove infinite values
    non.na.i = which(!is.na(z) & z != 0)
    z.non.na = z[non.na.i]
    
    ## Duplication
    z.dup = z.non.na[z.non.na > 0]
    z.dup = sample(c(-1, 1), length(z.dup), replace = TRUE) * z.dup
    if(!is.null(quant.int)){
      sd.df = plyr::ldply(quant.int, function(qi) {
        data.frame(quant = qi, sd.est = fdrtool::censored.fit(z.dup, quantile(abs(z.dup), 
                                 probs = qi, na.rm = TRUE))[5])
      })
      sd.e = localMax(sd.df$sd.est)$lM[1]
    } else {
      if(is.null(z.th)){stop("Either quant.int or z.th parameters are necessary.")}
      sd.e = fdrtool::censored.fit(z.dup, z.th[1])[5]
    }
    if (!is.null(ref.dist.weight)) {
        sd.df = with(sd.df, dplyr::arrange(sd.df, abs(sd.est - sd.e)))
        sd.e = sd.e + ref.dist.weight * 10 * sd(sd.df$sd.est[1:20])
    }
    res$sigma.est.dup = sd.e
    ## Deletion
    z.del = z.non.na[z.non.na < 0]
    z.del = sample(c(-1, 1), length(z.del), replace = TRUE) * z.del
    if(!is.null(quant.int)){
      sd.df = plyr::ldply(quant.int, function(qi) {
        data.frame(quant = qi, sd.est = fdrtool::censored.fit(z.del, quantile(abs(z.del), 
                                 probs = qi, na.rm = TRUE))[5])
      })
      sd.e = localMax(sd.df$sd.est)$lM[1]
    } else {
      if(length(z.th)>1){
        sd.e = fdrtool::censored.fit(z.del, z.th[2])[5]
      } else {
        sd.e = fdrtool::censored.fit(z.del, z.th[1])[5]
      }
    }
    if (!is.null(ref.dist.weight)) {
        sd.df = with(sd.df, dplyr::arrange(sd.df, abs(sd.est - sd.e)))
        sd.e = sd.e + ref.dist.weight * 10 * sd(sd.df$sd.est[1:20])
    }
    res$sigma.est.del = sd.e
    
    res$pval[non.na.i[z.non.na > 0]] = 2 * pnorm(-abs(z.dup), 0, res$sigma.est.dup)
    res$pval[non.na.i[z.non.na < 0]] = 2 * pnorm(-abs(z.del), 0, res$sigma.est.del)
    if (any(res$pval == 0, na.rm = TRUE)) 
        res$pval[which(res$pval == 0)] = .Machine$double.xmin
    res$qval = p.adjust(res$pval, method = "fdr")
    
    if (plot & any(!is.na(res$pval))) {
        pv = qv = NULL  ## Uglily appease R checks
        plot.df = data.frame(z = z, pv = res$pval, qv = res$qval)
        print(ggplot2::ggplot(plot.df[which(abs(plot.df$z) < quantile(abs(plot.df$z), 
            probs = 0.95, na.rm=TRUE) + 1), ], ggplot2::aes(x = z)) + ggplot2::geom_histogram() + 
            ggplot2::xlab("Z-score") + ggplot2::ylab("number of bins") + ggplot2::theme_bw())
        print(ggplot2::ggplot(plot.df, ggplot2::aes(x = pv)) + ggplot2::geom_histogram() + 
            ggplot2::xlab("P-value") + ggplot2::xlim(0, 1) + ggplot2::ylab("number of bins") + 
            ggplot2::theme_bw())
        print(ggplot2::ggplot(plot.df[which(plot.df$qv < 0.1), ], ggplot2::aes(x = cut(qv, 
            breaks = c(-Inf, 0.001, 0.01, 0.5, 0.1)))) + ggplot2::geom_bar() + ggplot2::xlab("Q-value") + 
            ggplot2::ylab("number of bins") + ggplot2::theme_bw())
    }
    
    return(res)
} 
