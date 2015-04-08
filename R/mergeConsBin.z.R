##' Merge consecutive abnormal bins into one larger event. Two consecutive bins
##' will be merged if the maximum (minimum) of their two Z-scores if significantly
##' lower (higher) than what is observed from a null distribution. The null distribution
##' is simulated by taking the maximum (minimum) values of two Normal values. The variance
##' of the simulated Normal is given as a parameter of the function.
##' @title Merge abnormal consecutive bins
##' @param res.df a data.frame with the Z-scores. Columns 'chr', 'start', 'end' and 'z'
##' are required.
##' @param fdr.th the False Discovery Rate threshold.
##' @param sd.null the estimated standard deviation of the Z-score null distribution.
##' Usually, computed during P-value estimation by 'fdrtool.quantile'.
##' @param nb.sim the number of simulated Z-scores for the P-value computation.
##' @return a data.frame similar to the input 'res.df' but with an extra 'nb.bin.cons'
##' column (the number of bin merged for each event).
##' @author Jean Monlong
##' @keywords internal
mergeConsBin.z <- function(res.df, fdr.th = 0.05, sd.null = 1, nb.sim = 1e+06) {

  col.mean = c("z", "fc", "mean.cov", "GCcontent", "lowC", "map")
  col.mean.log = c("pv", "qv")
  
  fun3 <- function(x, FUN = mean, log.x=FALSE) {
    if(log.x){x = log(x)}
    if (length(x) > 2) {
      res.x = FUN(x[2:(length(x)-1)])
    } else if (length(x) == 2) {
      res.x = FUN(x)
    } else {
      res.x = x
    }
    if(log.x){res.x = exp(res.x)}
    res.x
  }

  ## Compute P-value from an empirical null distribution
  compute.pv <- function(z, z.null, alt.greater = TRUE) {
    if (length(z) > 1) {
      z.f = factor(z)
      zn.c = cut(z.null, c(-Inf, levels(z.f), Inf), right = FALSE)
      zn.sum = table(zn.c)
      zn.cs = cumsum(zn.sum)[-length(zn.sum)]
      names(zn.cs) = levels(z.f)
      zn.cs.f = zn.cs[z.f]  ## Pb ? as.character, luckily not
      names(zn.cs.f) = names(z)
      if (alt.greater) {
        return(1 - (zn.cs.f/(length(z.null) + 1)))
      } else {
        return((zn.cs.f + 1)/(length(z.null) + 1))
      }
    } else {
      if (alt.greater) {
        pv = (sum(z <= z.null) + 1)/(length(z.null) + 1)
      } else {
        pv = (sum(z >= z.null) + 1)/(length(z.null) + 1)
      }
      names(pv) = names(z)
      return(pv)
    }
  }
  ## Simulate the empirical null distribution
  ## z.null = apply(rbind(rnorm(nb.sim, 0, sd.null), rnorm(nb.sim, 0, sd.null)), 2, sort)
  z.null = apply(rbind(sample(res.df$z, nb.sim, replace=TRUE), sample(res.df$z, nb.sim, replace=TRUE)), 2, sort)
  ## Compute P-values
  res.df = with(res.df, dplyr::arrange(res.df, chr, start))
  pvLink.f <- function(df) {
    z.link = apply(rbind(df$z[-1], df$z[-nrow(df)]), 2, sort)
    data.frame(pv.dup = compute.pv(z.link[1, ], z.null = z.null[1, ]), pv.del = compute.pv(z.link[2, ], z.null = z.null[2, ], alt.greater = FALSE))
  }
  link.df = lapply(unique(res.df$chr), function(chr.i)data.frame(chr=chr.i, pvLink.f(res.df[which(res.df$chr==chr.i),])))
  link.df = plyr::ldply(link.df, identity)
  ## Multiple test correction
  link.df$qv.dup = fdrtool::fdrtool(link.df$pv.dup, statistic = "pvalue", plot = FALSE, verbose = FALSE)$qval
  link.df$qv.del = fdrtool::fdrtool(link.df$pv.del, statistic = "pvalue", plot = FALSE, verbose = FALSE)$qval
  ## Annotate and merge bins
  link.annotate.f <- function(df) {
    ##cat(df$chr[1], "\n")
    link.df = link.df[which(link.df$chr == df$chr[1]), ]
    link.v = rep("none", nrow(df))
    link.v[c(FALSE, link.df$qv.dup <= fdr.th) | c(link.df$qv.dup <= fdr.th, FALSE) | (df$qv <= fdr.th & df$z > 0)] = "dup"
    link.v[c(FALSE, link.df$qv.del <= fdr.th) | c(link.df$qv.del <= fdr.th, FALSE) | (df$qv <= fdr.th & df$z < 0)] = "del"
    rl = rle(link.v)
    rlv = rl$values
    rlv[rlv != "none"] = paste(link.df$chr[1], rlv[rlv != "none"], 1:sum(rlv != "none"))
    df$link = rep(rlv, rl$lengths)
    return(df)
  }
  ##res.df = dplyr::do(dplyr::group_by(res.df, chr), link.annotate.f(.))
  res.df = lapply(unique(res.df$chr), function(chr.i)link.annotate.f(res.df[which(res.df$chr==chr.i),]))
  res.df = plyr::ldply(res.df, identity)
  res.df = res.df[which(res.df$link != "none"), ]
  if (nrow(res.df) > 0) {
    merge.event.f <- function(df.f) {
      df.o = with(df.f, data.frame(chr = head(chr,1), start = min(start), end = max(end), nb.bin.cons = nrow(df.f)))
      cbind(df.o,
            t(apply(df.f[, intersect(colnames(df.f), col.mean), drop = FALSE], 2, fun3)),
            t(apply(df.f[, intersect(colnames(df.f), col.mean.log), drop = FALSE], 2, fun3, log.x=TRUE))
            )
    }
    ##res.df = as.data.frame(dplyr::do(dplyr::group_by(res.df, link, chr), merge.event.f(.)))
    res.df = lapply(unique(res.df$link), function(l.i) merge.event.f(res.df[which(res.df$link==l.i),]))
    res.df = plyr::ldply(res.df, identity)
  } else {
    res.df = NULL
  }
  ##res.df$link = NULL
  return(res.df)
} 
