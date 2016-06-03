##' Detect abnormal bin from the Z-score distribution. A normal distribution is first
##' fitted to the Z-score distribution. P-values are computed from this estimated
##' null distribution and corrected for multiple testing. Eventually consecutive bins
##' with abnormal read counts can be merged.
##'
##' Two approaches can be used to define if a bin has abnormal threshold. By default ('sdest'), the null Normal distribution standard deviation is estimated by sequencially trimming the Z-score distribution and using an estimator for censored values. Once the Z-scores corresponding to the abnormal bins are trimmed out, the estimator reaches a plateau which is used as estimator for the null standard deviation. Using this parameter, P-values and Q-values are computed; abnormal bins are then defined by a user-defined FDR threshold on the Q-values. An alternative approach, 'consbins', looks at the distribution of consecutive bins to define the best threshold on the Z-scores. A wide range of thresholds are eplored. For each threshold, selected bins are stitched together if directly consecutive and the proportion of single  and pair bins is computed. With a loose value-many selected bins-, pairs of consecutive bins happen by chance. More stringent values decreases the proportion of pairs and increases the number of single bins until it reaches true calls that are more likely to be consecutive. The Z-score threshold is defined as the changepoint between random and true calls distribution. Eventually another version of 'sdest' is implemented but this time fitting two Gaussian distribution (centered in 0). This approach, 'sdest2N', is more suited when we suspect that the sample tested is not completely comparable to the reference samples. With the two Gaussian distribution a longer tail can be integrated in the null distribution, reducing the potential false calls in presence of a long-tail.
##'
##' Two approaches are available to merge bins with abnormal read coverage. 'stitch' simply stitches bins passing a user-defined significance threshold. In this approach, the stitching distance specifies the maximum distance between two bins that will be merged. By default the bin size is used, i.e. two abnormal bins will be merged if separated by maximum one bin. 'zscores' approach looks at the Z-score of two consecutive bins: if the minimum(maximum) is significantly higher(lower) than a simulated null distribution, these two bins will be merged to create a larger duplication(deletion).
##'
##' For cancer samples, 'min.normal.prop' can be reduced, e.g. to 0.6. Aneuploid can also be removed with 'aneu.chrs'. Function 'aneuploidy.flag' can help flagging aneuploid chromosomes.
##' @title Call abnormal bins
##' @param files.df a data.frame with the paths to different sample files (bin count, Z-scores, ..). Here columns 'z' and 'fc' are used to retrieve Z-scores and fold changes.
##' @param samp the name of the sample to analyze.
##' @param out.pdf the name of the output pdf file.
##' @param FDR.th the False Discovery Rate to use for the calls. 
##' @param merge.cons.bins how the bins should be merged. Default is 'stitch'. 'zscores' is another approch (see Details), 'no' means no bin merging.
##' @param stitch.dist the maximal distance between two calls to be merged into one (if 'merge.cons.bins="stitch"'). If NULL (default), the bin size + 1 is used.
##' @param z.th how the threshold for abnormal Z-score is chosen. Default is 'sdest' which will use 'FDR.th=' parameter as well. 'consbins' looks at the number of consecutive bins, see Details.
##' @param norm.stats the name of the file with the normalization statistics ('norm.stats' in 'tn.norm' function) or directly a 'norm.stats' data.frame.
##' @param min.normal.prop the minimum proportion of the regions expected to be normal. Default is 0.9. For cancers with many large aberrations, this number can be lowered. Maximum value accepted is 0.98 .
##' @param aneu.chrs the names of the chromosomes to remove because flagged as aneuploid. If NULL (default) all chromosomes are analyzed.
##' @param gc.df a data.frame with the GC content in each bin, for the Z-score normalization. Columns required: chr, start, end, GCcontent. If NULL (default), no normalization is performed.
##' @return a data.frame with columns
##' \item{chr, start, end}{the genomic region definition.}
##' \item{z}{the Z-score.}
##' \item{pv, qv}{the P-value and Q-value(~FDR).}
##' \item{fc}{the copy number estimate (if 'fc' was not NULL).}
##' \item{nb.bin.cons}{the number of consecutive bins (if the bins were merged, i.e. ' 'merge.cons.bins!='no'').}
##' \item{cn2.dev}{Copy number deviation from the reference.}
##' @author Jean Monlong
##' @export
call.abnormal.cov <- function(files.df=NULL, samp, out.pdf = NULL, FDR.th = 0.05, merge.cons.bins = c("stitch", "zscores", "cbs", "no"), stitch.dist=NULL, z.th = c("sdest", "consbins", "sdest2N"), norm.stats = NULL, min.normal.prop = 0.9, aneu.chrs = NULL, gc.df=NULL) {

  z.f = files.df$z[which(files.df$sample==samp)]
  if(!file.exists(z.f)) z.f = paste0(z.f, ".bgz")
  res.df = utils::read.table(z.f, header=TRUE, as.is=TRUE)
  fc.f = files.df$fc[which(files.df$sample==samp)]
  if(!file.exists(fc.f)) fc.f = paste0(fc.f, ".bgz")
  fc = utils::read.table(fc.f, header=TRUE, as.is=TRUE)
  ## Check the order is consistent in both files
  rand.ii = sample.int(nrow(res.df),100)
  if(nrow(res.df)!=nrow(fc) | any(res.df$start[rand.ii]!=fc$start[rand.ii]) | any(res.df$chr[rand.ii]!=fc$chr[rand.ii])){
    stop("Z-score and fold-change (fc) files are not in the same order. Maybe recompute them ?")
  }
  ##
  res.df$fc = fc$fc
  rm(fc)

  res.df = res.df[which(!is.na(res.df$z)),]
  
  if (!is.null(norm.stats)) {
    if (is.character(norm.stats) & length(norm.stats) == 1) {
      headers = utils::read.table(norm.stats, nrows = 1, as.is = TRUE)
      colC = rep("NULL", length(headers))
      colC[1:4] = c("character","integer","integer","numeric")
      norm.stats = utils::read.table(norm.stats, header = TRUE, colClasses = colC)
    }
    colnames(norm.stats)[4] = "mean.cov"
    res.df = merge(res.df, norm.stats, all.x=TRUE)
  }
  if(!is.null(gc.df)){
    res.df$z = z.norm(res.df, gc.df)
  }

  ## Remove aneuploid chromosomes
  if (!is.null(aneu.chrs)) {
    res.df = res.df[which(!(res.df$chr %in% aneu.chrs)), ]
  }

  if (!is.null(out.pdf)) {
    grDevices::pdf(out.pdf, 13, 10)
  }

  res.df = res.df[which(!is.na(res.df$z) & !is.infinite(res.df$z) & res.df$z!=0),]
  bin.width = stats::median(round(res.df$end - res.df$start + 1), na.rm=TRUE)
  ## Pvalue/Qvalue estimation
  if (all(is.na(res.df$z)))
    return(NULL)
  if (z.th[1] == "sdest") {
    if (min.normal.prop > 0.98) {
      stop("Maximum value accepted for 'min.normal.prop' is 0.98.")
    }
    fdr = fdrtool.quantile(res.df$z, quant.int = seq(min.normal.prop, 0.99, 0.01), plot = !is.null(out.pdf))
    res.df$pv = fdr$pval
    res.df$qv = fdr$qval

    ## Remove large aberrations
    ## aber.large = mergeConsBin.reduce(res.df[which(res.df$qv < 0.05), ], stitch.dist = 10 * bin.width)
    ## aber.large = subset(aber.large, end - start > 1e+07)
    ## if (nrow(aber.large) > 0) {
    ##   aber.gr = with(aber.large, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    ##   res.gr = with(res.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    ##   res.aber.large = res.df[which(GenomicRanges::overlapsAny(res.gr, aber.gr)), ]
    ##   res.df = res.df[which(!GenomicRanges::overlapsAny(res.gr, aber.gr)), ]
    ##   fdr = fdrtool.quantile(res.df$z, quant.int = seq(min.normal.prop, 0.99, 0.01), plot =  !is.null(out.pdf))
    ##   res.df$pv = fdr$pval
    ##   res.df$qv = fdr$qval
    ##   if (nrow(aber.large) > 0) {
    ##     res.df = rbind(res.df, res.aber.large)
    ##   }
    ## }

  } else if (z.th[1] == "consbins") {
    res.df = z.thres.cons.bins(res.df, plot = !is.null(out.pdf), pvalues = TRUE)$z.df
  } else if (z.th[1] == "sdest2N") {
    fdr = fdrtool.quantile.2N(res.df$z, plot = !is.null(out.pdf))
    res.df$pv = fdr$pval
    res.df$qv = fdr$qval
  } else {
    stop("'z.th=': available thresholding approaches are : 'stitch', 'zscores'.")
  }

  if (merge.cons.bins[1] != "no") {
    if(is.null(stitch.dist)) {
      stitch.dist = bin.width + 1
    }

    if (merge.cons.bins[1] == "stitch") {
      res.df = res.df[which(res.df$qv < FDR.th), ]
      res.df = mergeConsBin.reduce(res.df, stitch.dist = stitch.dist)
    } else if (merge.cons.bins[1] == "zscores") {
      res.df = mergeConsBin.z(res.df, fdr.th = FDR.th, sd.null = max(c(fdr$sigma.est.dup,fdr$sigma.est.del)), stitch.dist = stitch.dist)
    } else if (merge.cons.bins[1] == "cbs") {
      res.df = mergeConsBin.cbs(res.df, pv.th = FDR.th)
    } else {
      stop("'merge.cons.bins=' : available bin merging approaches are : 'stitch', 'zscores'.")
    }

    if (nrow(res.df) > 0 & !is.null(out.pdf)) {
      nb.bin.cons = NULL  ## Uglily appease R checks
      print(ggplot2::ggplot(res.df, ggplot2::aes(x = factor(nb.bin.cons))) +
              ggplot2::geom_bar() + ggplot2::ylab("number of bins") + ggplot2::xlab("number of consecutive abnormal bins") +
                ggplot2::theme_bw())

      if (any(colnames(res.df) == "fc") & sum(res.df$nb.bin.cons > 2 & res.df$fc < 2.5) > 3) {
        print(ggplot2::ggplot(subset(res.df, nb.bin.cons > 2), ggplot2::aes(x = 2 * fc)) + ggplot2::geom_histogram() + ggplot2::theme_bw() + ggplot2::ylab("number of bins") + ggplot2::xlab("copy number estimate") + ggplot2::ggtitle("At least 3 consecutive abnormal bins") + ggplot2::xlim(0, 5))
      }
    }
  }

  if (!is.null(out.pdf)) {
    grDevices::dev.off()
  }

  ## Deviation from copy number 2
  res.df$cn2.dev = 2*abs(res.df$fc - 1)
  res.df$cn = round(2*res.df$fc)
  if(nrow(res.df)>0){
    res.df$prop.single.bin = mean(res.df$nb.bin.cons==1)
  } else {
    res.df$prop.single.bin = numeric(0)
  }

  if (nrow(res.df) > 0 & merge.cons.bins[1] != "no") {
    return(data.frame(sample = samp, res.df, stringsAsFactors = FALSE))
  } else if (any(res.df$qv <= FDR.th, na.rm = TRUE)) {
    return(data.frame(sample = samp, res.df[which(res.df$qv <= FDR.th), ], stringsAsFactors = FALSE))
  } else {
    return(NULL)
  }
}
