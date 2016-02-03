##' Compute P-values and Q-values from the Z-score distribution. Here the Z-score distribution is modeled as a Normal centered in 0. The variance is fitted at different level of tail-trimming.
##' @title P-values estimation
##' @param z a vector with the Z-scores
##' @param quant.int a vector with the quantiles to test for the estimation of the
##' null normal variance.
##' @param plot should some graphs be displayed. Default if FALSE.
##' @param z.th the Z-score threshold to use to estimate the null normal variance, if 'quant.int=NULL'. If one value the same threshold is used for both duplication and deletion; if a vector, the threshold for duplication / deletions.
##' @return a list with
##' \item{pval}{the vector of P-values}
##' \item{qval}{the vector of Q-values / FDR estimates}
##' \item{sigma.est.dup}{the estimated null distribution variance for positive Z-scores}
##' \item{sigma.est.del}{the estimated null distribution variance for negative Z-scores}
##' @author Jean Monlong
##' @keywords internal
fdrtool.quantile <- function(z, quant.int = seq(0.4, 1, 0.02), plot = TRUE, z.th=NULL) {
  localMax <- function(x, min.max.prop = 0.1) {
    ## Remove extreme outliers
    if(any(x > 10 * sort(x, decreasing=TRUE)[2])){
      x = x[which(x < 10 * sort(x, decreasing=TRUE)[2])]
    }
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
  res$sigma.est.del = sd.e

  res$pval[non.na.i[z.non.na > 0]] = 2 * pnorm(-abs(z.dup), 0, res$sigma.est.dup)
  res$pval[non.na.i[z.non.na < 0]] = 2 * pnorm(-abs(z.del), 0, res$sigma.est.del)
  if (any(res$pval == 0, na.rm = TRUE))
  res$pval[which(res$pval == 0)] = .Machine$double.xmin
  res$qval = p.adjust(res$pval, method = "fdr")

  if (plot & any(!is.na(res$pval))) {
    pv = qv = ..density.. = y = NULL  ## Uglily appease R checks
    plot.df = data.frame(z = z, pv = res$pval, qv = res$qval)

    z.lim = c(-res$sigma.est.del, res$sigma.est.dup)*ifelse(mean(res$pval<.01)>.1,8,5)
    null.df = data.frame(y=c(dnorm(seq(z.lim[1],0,.05),0,res$sigma.est.del),dnorm(seq(0,z.lim[2],.05),0,res$sigma.est.dup)), z=c(seq(z.lim[1],0,.05),seq(0,z.lim[2],.05)))
    null.df$y = null.df$y * mean(z> -4*res$sigma.est.del & z<4*res$sigma.est.dup)

    print(ggplot2::ggplot(plot.df, ggplot2::aes(x = z)) +
          ggplot2::geom_histogram(ggplot2::aes(y=..density..)) +
          ggplot2::xlab("Z-score") + ggplot2::ylab("number of bins") + ggplot2::theme_bw() +
          ggplot2::geom_line(ggplot2::aes(y=y), data=null.df, linetype=2, colour="red") +
          ggplot2::xlim(z.lim))

    print(ggplot2::ggplot(plot.df, ggplot2::aes(x = pv, fill=cut(qv, breaks = c(-Inf, 0.001, 0.01, 0.5, 0.1,1)))) + ggplot2::geom_histogram() +
          ggplot2::xlab("P-value") + ggplot2::xlim(0, 1) + ggplot2::ylab("number of bins") +
          ggplot2::scale_fill_hue(name="Q-value") +
          ggplot2::theme_bw() + ggplot2::theme(legend.position="bottom"))
  }

  return(res)
}
