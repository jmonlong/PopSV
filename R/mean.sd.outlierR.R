##' Computes the mean and standard deviation after removing outlier groups.
##' \code{outliers} package and Grubbs approach are used.
##' @title Mean and standard diviation after outlier removal
##' @param x a vector of numeric values.
##' @param pv.max.ol the maximum P-value for the Grubbs test. Any outlier
##' with a P-value lower will be considered as outlier.
##' @return a vector similar to 'x' but with NAs for outliers.
##' @author Jean Monlong
##' @keywords internal
mean.sd.outlierR <- function(x, pv.max.ol = 1e-06) {
  
  sd.mad <- function(x) {
    if (all(x == 0, na.rm = TRUE)) {
      return(sd(rpois(length(x), 1)))
    }
    sd.res = mad(x, na.rm = TRUE)
    if (sd.res == 0) {
      return(sd(rpois(length(x), 1)))
    } else {
      return(sd.res)
    }
  }
  trim.mean <- function(x, probs=c(.4,.6)){
    qq = quantile(x,probs=probs, na.rm=TRUE)
    x[x<qq[1] | x>qq[2]] = NA
    m = mean(x, na.rm=TRUE)
    if(m==0){
      return(1)
    }
    return(m)
  }
  fit2norm.msd <- function(z) {
    mix.obj <- function(p, x) {
      e <- p[1] * dnorm((x-p[2])/p[4])/p[4] + (1 - p[1]) * dnorm((x-p[3])/p[5])/ p[5]
      if (any(e <= 0, na.rm = TRUE) | p[1] < 0 | p[1] > 1) 
      Inf else -sum(log(e))
    }
    lmix2a <- deriv(~-log(p * dnorm((x-m1)/s1)/s1 + (1 - p) * dnorm((x-m2)/s2)/s2), c("p", "m1","m2", "s1", "s2"), function(x, p, m1, m2, s1, s2) NULL)
    mix.gr <- function(pa, x) {
      p <- pa[1]
      m1 <- pa[2]
      m2 <- pa[3]
      s1 <- pa[4]
      s2 <- pa[5]
      colSums(attr(lmix2a(x, p, m1, m2, s1, s2), "gradient"))
    }
    p0 = c(p=.9, m1=median(z), m2=mean(z), s1=mad(z)/2, s2=sd(z))
    results = optim(p0, mix.obj, mix.gr, x = z)
    if (results$par[1] < 0.5) {
      results$par[1] = 1 - results$par[1]
      results$par[2:3] = results$par[3:2]
      results$par[4:5] = results$par[5:4]
    }
    results
  }
  grubbs.sw <- function(x, type = 10, opposite = FALSE, two.sided = TRUE, max.pv = 0.01, max.step = 10) {
    grubbs.t <- function(x, type = 10, opposite = FALSE, two.sided = TRUE) {
      gt = outliers::grubbs.test(x, type, opposite, two.sided)
      gt$outlier = outliers::outlier(x)
      gt$outlier.i = which(outliers::outlier(x, logical = TRUE))[1]
      return(gt)
    }
    pv = outliers = NULL
    step = 1
    continue = TRUE
    while (continue & step <= max.step) {
      gt = grubbs.t(as.numeric(x), type, opposite, two.sided)
      if (gt$p.value > max.pv | all(is.na(gt$statistic))) {
        continue = FALSE
      } else {
        pv = c(pv, gt$p.value)
        outliers = c(outliers, gt$outlier.i)
        x[outliers[step]] = NA
      }
      step = step + 1
    }
    return(list(x = x, pv = pv, outliers = outliers, max.pv = max.pv))
  }

  if(all(is.na(x))){
    return(list(m=NA,sd=NA, nb.remove=NA))
  }
  
  xs = sort(x)
  dd = diff(xs)
  gt = grubbs.sw(dd, max.pv = pv.max.ol)
  if (!is.null(gt$outliers)) {
    cut.pos = c(0, sort(gt$outliers), length(xs))
    gp.sizes = diff(cut.pos)
    cm = which.max(gp.sizes)
    x.rm = rep(NA, length(xs))
    x.rm[(cut.pos[cm] + 1):cut.pos[cm + 1]] = xs[(cut.pos[cm] + 1):cut.pos[cm + 1]]
  } else {
    x.rm = xs
  }
  return(list(m = max(1,trim.mean(x.rm)), sd = max(1,sd.mad(x.rm)), nb.remove = length(x) - sum(!is.na(x.rm))))
  ## Mixture of Gaussian, more flexible/robust but 20 times slower...
  ##msd = fit2norm.msd(na.omit(x.rm))
  ##return(list(m = max(1,msd$par[2]), sd = max(1,msd$par[4]), nb.remove = length(x) - sum(!is.na(x.rm))))
} 
