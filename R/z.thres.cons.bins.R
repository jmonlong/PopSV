##' Find a threshold on the Z-scores associated with how many consecutive bins are found.
##' @title Z-score thresholding using bin consecutiveness.
##' @param z.df a data.frame with at least 'chr', 'start', 'end' and 'z' columns.
##' @return a list with:
##' \item{dup.th}{the threshold for positive Z-scores (duplication signal)}
##' \item{del.th}{the threshold for negative Z-scores (deletion signal)}
##' \item{sigma.est.dup}{the estimated null distribution variance for positive Z-scores}
##' \item{sigma.est.del}{the estimated null distribution variance for negative Z-scores}
##' @author Jean Monlong
##' @keywords internal
z.thres.cons.bins <- function(z.df) {

  bin.w = round(stats::median(z.df$end - z.df$start + 1, na.rm = TRUE))
  cons.dist.f <- function(df) {
    if(nrow(df)==0){return(data.frame(nbc = NA, n = NA, p = NA))}
    gr = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    gr.r = GenomicRanges::reduce(gr, min.gapwidth = 2)
    nbc.t = table(round(GenomicRanges::width(gr.r)/bin.w))
    data.frame(nbc = as.numeric(names(nbc.t)), n = as.numeric(nbc.t), p = as.numeric(nbc.t)/sum(nbc.t))
  }
  ## cumLocalMax <- function(x, y, min = FALSE) {
  ##   rec.max <- function(xx) {
  ##     if (length(xx) > 1) {
  ##       xx.p = rec.max(xx[-length(xx)])
  ##       if (min) {
  ##         return(c(xx.p, min(c(utils::tail(xx.p, 1), xx[length(xx)]))))
  ##       } else {
  ##         return(c(xx.p, max(c(utils::tail(xx.p, 1), xx[length(xx)]))))
  ##       }
  ##     } else {
  ##       return(xx)
  ##     }
  ##   }
  ##   cmax = rec.max(y)
  ##   cmax.rle = rle(diff(cmax))
  ##   if (all(cmax.rle$values != 0)) {
  ##     return(list(x = utils::tail(x, 1), y = utils::tail(y, 1)))
  ##   }
  ##   rle.i = which(cmax.rle$values == 0)[which.max(cmax.rle$lengths[cmax.rle$values ==
  ##                   0])]
  ##   x.i = sum(cmax.rle$lengths[1:rle.i]) - cmax.rle$lengths[rle.i] + 1
  ##   return(list(x = x[x.i], y = y[x.i]))
  ## }
  localMax <- function(x, y, min.max.prop = 0.1, loc.max = TRUE) {
    d = data.frame(x = x, y = y)
    d = dplyr::arrange(d, x)
    my = max(d$y)
    if (!loc.max) {
      d$y = my - d$y
      my = max(d$y)
    }
    im = 1 + which(diff(sign(diff(d$y))) == -2)
    max.id = im[which(d$y[im] >= min.max.prop * my)]
    if (length(max.id) == 0) {
      return(list(lM = NA, h = NA))
    }
    max.id.o = max.id[order(d$y[max.id], decreasing = TRUE)]
    return(list(lM = d$x[max.id.o], h = d$y[max.id.o]/my))
  }
  find.th <- function(df, z.int = seq(1, 20, 0.2)) {
    nbcc.df = do.call(rbind, lapply(z.int, function(z.th) {
      df.th = df[which(abs(df$z) > z.th), ]
      df.th = cons.dist.f(df.th)
      df.th$z.th = z.th
      return(df.th)
    }))
    ## ggplot(subset(nbcc.df, nbc<3), aes(x=z.th, y=p)) + geom_line() + facet_grid(nbc~., scales='free')
    z.th = suppressWarnings(c(min(localMax(nbcc.df$z.th[which(nbcc.df$nbc==1)], nbcc.df$p[which(nbcc.df$nbc==1)], loc.max=FALSE)$lM), sapply(2, function(nbc.i)min(localMax(nbcc.df$z.th[which(nbcc.df$nbc==nbc.i)], nbcc.df$p[which(nbcc.df$nbc==nbc.i)])$lM))))
    ## z.th = c(cumLocalMax(nbcc.df$z.th[which(nbcc.df$nbc == 1)], nbcc.df$p[which(nbcc.df$nbc == 1)], min = TRUE)$x, cumLocalMax(nbcc.df$z.th[which(nbcc.df$nbc == 2)], nbcc.df$p[which(nbcc.df$nbc == 2)])$x)
    if(all(is.na(z.th))) return(max(abs(df$z)))
    mean(z.th, na.rm = TRUE)
  }

  ## Split between duplication/deletion signal
  dup.df = z.df[which(z.df$z > 0), ]
  del.df = z.df[which(z.df$z < 0), ]

  ## Find threshold; second run scan with more resolution.
  dup.th = find.th(dup.df, z.int = seq(2, stats::quantile(dup.df$z, probs = 0.999) * 2, 0.2))
  ## dup.th = find.th(dup.df, seq(dup.th-.5, dup.th+.5, .01))
  ## dup.th = find.th(dup.df, seq(dup.th-.1, dup.th+.1, .005))
  del.th = find.th(del.df, z.int = seq(2, stats::quantile(abs(del.df$z), probs = 0.999) * 2, 0.2))
  ## del.th = find.th(del.df, seq(del.th-.5, del.th+.5, .01))
  ## del.th = find.th(del.df, seq(del.th-.1, del.th+.1, .005))

  ## P-value computation
  fdr = fdrtool.quantile(z.df$z, quant.int = NULL, z.th=c(dup.th, del.th))
  
  return(list(dup.th = dup.th, del.th = del.th, sigma.est.dup=fdr$sigma.est.dup, sigma.est.del=fdr$sigma.est.del))
}
