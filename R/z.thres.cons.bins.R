##' Find a threshold on the Z-scores associated with how many consecutive bins are found.
##' @title Z-score thresholding using bin consecutiveness.
##' @param z.df a data.frame with at least 'chr', 'start', 'end' and 'z' columns.
##' @param plot should some graphs be displayed. Default if FALSE.
##' @param pvalues should the P-values be deduced from these thresholds. Default is FALSE. Note: the P-values and Q-values derive from custom 'recipes' and should be taken as informative only, for further filtering or priorization.
##' @return a list with:
##' \item{dup.th}{the threshold for positive Z-scores (duplication signal)}
##' \item{del.th}{the threshold for negative Z-scores (deletion signal)}
##' \item{nb.ab.bins}{the number of abnormal bins}
##' \item{prop.ab.bins}{the proportion of abnormal bins}
##' \item{z.df}{a data.frame, similar to the input 'z.df' with new columns: 'abnormal' for bins with abnormal Z-scores and potentially 'pv'/'qv' for P-values/Q-values.}
##' @author Jean Monlong
##' @keywords internal
z.thres.cons.bins <- function(z.df, plot = FALSE, pvalues = FALSE) {

  bin.w = round(stats::median(z.df$end - z.df$start + 1, na.rm = TRUE))
  cons.dist.f <- function(df) {
    if(nrow(df)==0){return(data.frame(nbc = NA, n = NA, p = NA))}
    gr = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    gr.r = GenomicRanges::reduce(gr, min.gapwidth = 2)
    nbc.t = table(round(GenomicRanges::width(gr.r)/bin.w))
    data.frame(nbc = as.numeric(names(nbc.t)), n = as.numeric(nbc.t), p = as.numeric(nbc.t)/sum(nbc.t))
  }
  cumLocalMax <- function(x, y, min = FALSE) {
    rec.max <- function(xx) {
      if (length(xx) > 1) {
        xx.p = rec.max(xx[-length(xx)])
        if (min) {
          return(c(xx.p, min(c(utils::tail(xx.p, 1), xx[length(xx)]))))
        } else {
          return(c(xx.p, max(c(utils::tail(xx.p, 1), xx[length(xx)]))))
        }
      } else {
        return(xx)
      }
    }
    cmax = rec.max(y)
    cmax.rle = rle(diff(cmax))
    if (all(cmax.rle$values != 0)) {
      return(list(x = utils::tail(x, 1), y = utils::tail(y, 1)))
    }
    rle.i = which(cmax.rle$values == 0)[which.max(cmax.rle$lengths[cmax.rle$values ==
                    0])]
    x.i = sum(cmax.rle$lengths[1:rle.i]) - cmax.rle$lengths[rle.i] + 1
    return(list(x = x[x.i], y = y[x.i]))
  }
  localMax <- function(x, y = NULL, min.max.prop = 0.1, loc.max = TRUE) {
    if (is.null(y)) {
      d = stats::density(x, na.rm = TRUE)
    } else {
      d = data.frame(x = x, y = y)
      d = dplyr::arrange(d, x)
    }
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

  ## Annotate data.frames
  dup.df$abnormal = dup.df$z > dup.th
  del.df$abnormal = del.df$z < -del.th
  z.df = rbind(dup.df, del.df)

  ## P-value computation
  if (pvalues) {
    fdr = fdrtool.quantile(z.df$z, quant.int = NULL, plot = FALSE, z.th=c(dup.th, del.th))
    z.df$pv = fdr$pval
    z.df$qv = fdr$qval
  }

  ## Some statistics on the calls
  nb.ab.bins = sum(z.df$abnormal)
  prop.ab.bins = mean(z.df$abnormal)

  ## Z-score distribution with thresholds
  if (plot) {
    z = pv = qv = abnormal = NULL  ## Uglily appease R checks
    print(ggplot2::ggplot(z.df, ggplot2::aes(x = z)) + ggplot2::geom_histogram() +
          ggplot2::xlim(-20, 20) + ggplot2::geom_vline(xintercept = c(-del.th, dup.th), linetype = 2) + ggplot2::theme_bw())
    if (pvalues) {
      print(ggplot2::ggplot(z.df, ggplot2::aes(x = pv, fill=abnormal)) + ggplot2::geom_histogram() +
             ggplot2::xlim(0, 1) + ggplot2::theme_bw())
    }
  }

  return(list(dup.th = dup.th, del.th = del.th, nb.ab.bins = nb.ab.bins, prop.ab.bins = prop.ab.bins, z.df = z.df[which(z.df$abnormal), ]))
}
