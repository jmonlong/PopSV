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
##' @param stitch.dist the stitching distance, i.e. the maximum distance at which two bins will be merged.
##' @return a data.frame similar to the input 'res.df' but with an extra 'nb.bin.cons'
##' column (the number of bin merged for each event).
##' @author Jean Monlong
##' @keywords internal
mergeConsBin.z <- function(res.df, fdr.th = 0.05, sd.null = 1, nb.sim = 1e+06, stitch.dist=1e4) {

  col.mean = c("z", "fc", "mean.cov", "GCcontent", "lowC", "map")
  col.min = c("pv", "qv")

  fun3 <- function(x, FUN3 = mean) {
    if (length(x) > 2) {
      res.x = FUN3(x[2:(length(x)-1)], na.rm=TRUE)
    } else if (length(x) == 2) {
      res.x = FUN3(x, na.rm=TRUE)
    } else {
      res.x = x
    }
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
  z.null = apply(rbind(stats::rnorm(nb.sim, 0, sd.null), stats::rnorm(nb.sim, 0, sd.null)), 2, sort)
  ## Compute P-values
  res.df = with(res.df, dplyr::arrange(res.df, chr, start))
  pvLink.f <- function(df) {
    z.link = apply(rbind(df$z[-1], df$z[-nrow(df)]), 2, sort)
    data.frame(start=df$start[-nrow(df)],end=df$end[-1],pv.dup = compute.pv(z.link[1, ], z.null = z.null[1, ]), pv.del = compute.pv(z.link[2,], z.null = z.null[2, ], alt.greater = FALSE))
  }
  link.df = lapply(unique(res.df$chr), function(chr.i)data.frame(chr=chr.i, pvLink.f(res.df[which(res.df$chr==chr.i),])))
  link.df = do.call(rbind, link.df)
  ## Multiple test correction
  link.df$qv.dup = fdrtool::fdrtool(link.df$pv.dup, statistic = "pvalue", plot = FALSE, verbose = FALSE)$qval
  link.df$qv.del = fdrtool::fdrtool(link.df$pv.del, statistic = "pvalue", plot = FALSE, verbose = FALSE)$qval
  ## Annotate and merge bins

  gr.f = with(res.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
  res.df$red.i = NA

  ## Merge duplications
  dup.gr = with(link.df[which(link.df$qv.dup <= fdr.th), ], GenomicRanges::GRanges(chr, IRanges::IRanges(start,end)))
  dup.gr = GenomicRanges::reduce(dup.gr, min.gapwidth = stitch.dist)
  dup.gr = with(res.df[which(res.df$qv <= fdr.th & res.df$z>0), ], c(dup.gr, GenomicRanges::GRanges(chr, IRanges::IRanges(start,end))))
  dup.gr = GenomicRanges::reduce(dup.gr, min.gapwidth = stitch.dist)
  ol.dup = GenomicRanges::findOverlaps(gr.f, dup.gr)
  res.df$red.i[S4Vectors::queryHits(ol.dup)] = paste0("dup", S4Vectors::subjectHits(ol.dup))
  ## Merge deletions
  del.gr = with(link.df[which(link.df$qv.del <= fdr.th), ], GenomicRanges::GRanges(chr, IRanges::IRanges(start,end)))
  del.gr = GenomicRanges::reduce(del.gr, min.gapwidth = stitch.dist)
  del.gr = with(res.df[which(res.df$qv <= fdr.th & res.df$z<0), ], c(del.gr, GenomicRanges::GRanges(chr, IRanges::IRanges(start,end))))
  del.gr = GenomicRanges::reduce(del.gr, min.gapwidth = stitch.dist)
  ol.del = GenomicRanges::findOverlaps(gr.f, del.gr)
  res.df$red.i[S4Vectors::queryHits(ol.del)] = paste0("del", S4Vectors::subjectHits(ol.del))

  merge.event.f <- function(df.f) {
    df.o = with(df.f, data.frame(start = min(start), end = max(end), nb.bin.cons = nrow(df.f)))
    cbind(df.o,
          t(apply(df.f[, intersect(colnames(df.f), col.mean), drop = FALSE], 2, fun3)),
          t(apply(df.f[, intersect(colnames(df.f), col.min), drop = FALSE], 2, fun3, FUN3=min))
          )
  }

  red.i = chr = . = NULL  ## Uglily appease R checks
  res.df = res.df[which(!is.na(res.df$red.i)),]
  res.df = as.data.frame(dplyr::do(dplyr::group_by(res.df, red.i, chr), merge.event.f(.)))
  res.df$red.i = NULL
  return(res.df)
}
