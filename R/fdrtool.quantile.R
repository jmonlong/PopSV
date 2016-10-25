##' Compute the null variances. The Z-score distribution is modeled as a Normal centered in 0. The variance is fitted at different level of tail-trimming.
##' @title P-values estimation
##' @param z a vector with the Z-scores
##' @param quant.int a vector with the quantiles to test for the estimation of the
##' null normal variance.
##' @param z.th the Z-score threshold to use to estimate the null normal variance, if 'quant.int=NULL'. If one value the same threshold is used for both duplication and deletion; if a vector, the threshold for duplication / deletions.
##' @return a list with
##' \item{sigma.est.dup}{the estimated null distribution variance for positive Z-scores}
##' \item{sigma.est.del}{the estimated null distribution variance for negative Z-scores}
##' @author Jean Monlong
##' @keywords internal
fdrtool.quantile <- function(z, quant.int = seq(0.4, 1, 0.02), z.th=NULL) {
  localMax <- function(x, min.max.prop = 0.1) {
    ## Remove extreme outliers
    x = x[which(x<100)]
    if(any(x > 20 * stats::median(x, na.rm=TRUE))){
      x = x[which(x < 20 * stats::median(x, na.rm=TRUE))]
    }
    d = stats::density(x, na.rm = TRUE)
    im = 1 + which(diff(sign(diff(d$y))) == -2)
    my = max(d$y[im])
    max.id = im[which(d$y[im] >= min.max.prop * my)]
    max.id.o = max.id[order(d$y[max.id], decreasing = TRUE)]
    return(list(lM = d$x[max.id.o], h = d$y[max.id.o]/my))
  }

  res = list(sigma.est.dup = NA, sigma.est.del = NA)
  z = z[which(!is.infinite(z) & !is.na(z) & z != 0)]

  ## Duplication
  z.dup = z[z > 0]
  z.dup = sample(c(-1, 1), length(z.dup), replace = TRUE) * z.dup
  if(!is.null(quant.int)){
    sd.df = plyr::ldply(quant.int, function(qi) {
      data.frame(quant = qi, sd.est = fdrtool::censored.fit(z.dup, stats::quantile(abs(z.dup),
                               probs = qi, na.rm = TRUE))[5])
    })
    sd.e = localMax(sd.df$sd.est)$lM[1]
  } else {
    if(is.null(z.th)){stop("Either quant.int or z.th parameters are necessary.")}
    sd.e = fdrtool::censored.fit(z.dup, z.th[1])[5]
  }
  res$sigma.est.dup = sd.e

  ## Deletion
  z.del = z[z < 0]
  z.del = sample(c(-1, 1), length(z.del), replace = TRUE) * z.del
  if(!is.null(quant.int)){
    sd.df = plyr::ldply(quant.int, function(qi) {
      data.frame(quant = qi, sd.est = fdrtool::censored.fit(z.del, stats::quantile(abs(z.del),
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

  return(res)
}
