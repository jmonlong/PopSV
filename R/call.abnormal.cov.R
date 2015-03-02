##' Detect abnormal bin from the Z-score distribution. A normal distribution is first
##' fitted to the Z-score distribution. P-values are computed from this estimated
##' null distribution and corrected for multiple testing. Eventually consecutive bins
##' with abnormal read counts can be merged.
##'
##' Two approaches can be used to define if a bin has abnormal threshold. By default ('sdest'), the null Normal distribution standard deviation is estimated by sequencially trimming the Z-score distribution and using an estimator for censored values. Once the Z-scores corresponding to the abnormal bins are trimmed out, the estimator reaches a plateau which is used as estimator for the null standard deviation. Using this parameter, P-values and Q-values are computed; abnormal bins are then defined by a user-defined FDR threshold on the Q-values. An alternative approach, 'consbins', looks at the distribution of consecutive bins to define the best threshold on the Z-scores. A wide range of thresholds are eplored. For each threshold, selected bins are stitched together if directly consecutive and the proportion of single  and pair bins is computed. With a loose value-many selected bins-, pairs of consecutive bins happen by chance. More stringent values decreases the proportion of pairs and increases the number of single bins until it reaches true calls that are more likely to be consecutive. The Z-score threshold is defined as the changepoint between random and true calls distribution.
##'
##' Two approaches are available to merge bins with abnormal read coverage. 'stitch' simply stitches bins passing a user-defined significance threshold. In this approach, the stitching distance specifies the maximum distance between two bins that will be merged. By default the bin size is used, i.e. two abnormal bins will be merged if separated by maximum one bin. 'zscores' approach looks at the Z-score of two consecutive bins: if the minimum(maximum) is significantly higher(lower) than a simulated null distribution, these two bins will be merged to create a larger duplication(deletion).
##'
##' DETAILS ON MIXING DIFFERENT SAMPLES / TUMORS
##' @title Call abnormal bins
##' @param z the name of the file with the Z-scores for all samples OR or a
##' data.frame with the Z-scores for all samples. 
##' @param files.df a data.frame with the paths to different sample files (bin count, Z-scores, ..). Here column 'z' and 'fc' are used to retrieve Z-scores and fold changes.
##' @param samp the name of the sample to analyze.
##' @param out.pdf the name of the output pdf file.
##' @param FDR.th the False Discovery Rate to use for the calls. Further filtering
##' can always be performed.
##' @param merge.cons.bins how the bins should be merged. Default is 'stitch'. 'zscores' is another approch (see Details), 'no' means no bin merging.
##' @param z.th how the threshold for abnormal Z-score is chosen. Default is 'sdest' which will use 'FDR.th=' parameter as well. 'consbins' looks at the number of consecutive bins, see Details.
##' @param fc the name of the file with the copy number estimates for all samples OR
##' or a data.frame with the copy number estimates for all samples. 
##' @param norm.stats the name of the file with the normalization statistics ('norm.stats' in 'tn.norm' function) or directly a 'norm.stats' data.frame.
##' @param d.max.max the maximum correlation of the last supporting bin. 
##' @param min.normal.prop the minimum proportion of the regions expected to be normal. Default is 0.5. For cancers with many large aberrations, this number can be lowered. Maximum value accepted is 0.98 . 
##' @param aneu.chrs the names of the chromosomes to remove because flagged as aneuploid. If NULL (default) all chromosomes are analyzed.
##' @param ref.dist.weight the weight (value between 0 and 1) based on the distance to the reference samples.
##' @return a data.frame with columns
##' \item{chr, start, end}{the genomic region definition}
##' \item{z}{the Z-score}
##' \item{pv, qv}{the P-value and Q-value(~FDR)}
##' \item{fc}{the copy number estimate (if 'fc' was not NULL).}
##' \item{nb.bin.cons}{the number of consecutive bins (if the bins were merged;
##' 'merge.cons.bins!='no'')}
##' \item{cn2.dev}{Copy number deviation from the reference }
##' @author Jean Monlong
##' @export
call.abnormal.cov <- function(z=NULL, files.df=NULL, samp, out.pdf = NULL, FDR.th = 0.05, merge.cons.bins = c("stitch", 
    "zscores", "no"), z.th = c("sdest", "consbins", "sdest2N"), fc = NULL, norm.stats = NULL, 
    d.max.max = 0.5, min.normal.prop = 0.9, aneu.chrs = NULL, ref.dist.weight = NULL) {

  if(!is.null(z)){
    ## load Z-scores and FC coefficients
    if (is.character(z) & length(z) == 1) {
      headers = read.table(z, nrows = 1, as.is = TRUE)
      colC = rep("NULL", length(headers))
      if (!all(c("chr", "start", "end", samp) %in% headers)) {
        stop("Columns missing in Z file. Check that 'chr', 'start', 'end' and the sample column are present.")
      }
      colC[headers %in% c("chr", "start", "end", samp)] = c("character", rep("integer", 
            2), "numeric")
      res.df = read.table(z, header = TRUE, colClasses = colC)
    } else {
      res.df = z[, c("chr", "start", "end", samp)]
      rm(z)
    }
    colnames(res.df)[4] = "z"
    if (!is.null(fc)) {
      if (is.character(fc) & length(fc) == 1) {
        headers = read.table(fc, nrows = 1, as.is = TRUE)
        colC = rep("NULL", length(headers))
        if (all(headers != samp)) {
          stop("Columns missing in FC file. Check that 'chr', 'start', 'end' and the sample column are present.")
        }
        colC[headers == samp] = "numeric"
        fc = read.table(fc, header = TRUE, colClasses = colC)
    }
      res.df$fc = fc[, make.names(samp)]
    rm(fc)
    }
  } else if(!is.null(files.df)) {
    z.f = subset(files.df, sample==samp)$z
    if(!file.exists(z.f)) z.f = paste0(z.f, ".bgz")
    res.df = read.table(z.f, as.is=TRUE, header=TRUE)
    fc.f = subset(files.df, sample==samp)$fc
    if(!file.exists(fc.f)) fc.f = paste0(fc.f, ".bgz")
    fc = read.table(fc.f, as.is=TRUE, header=TRUE)
    res.df$fc = fc$fc
    rm(fc)
  } else {
    stop("Either 'z' or 'files.df' pameter must be assigned.")
  }
  
  if (!is.null(norm.stats)) {
        if (is.character(norm.stats) & length(norm.stats) == 1) {
            headers = read.table(norm.stats, nrows = 1, as.is = TRUE)
            colC = ifelse(headers == "d.max", "numeric", "NULL")
            d.max = read.table(norm.stats, header = TRUE, colClasses = colC)
        }
        res.df$d.max = d.max$d.max
        rm(d.max)
        res.df = subset(res.df, !is.na(d.max) & d.max != 1 & d.max < d.max.max)
    }
    
    ## Remove aneuploid chromosomes
    if (!is.null(aneu.chrs)) {
        res.df = res.df[which(!(res.df$chr %in% aneu.chrs)), ]
    }
    
    if (!is.null(out.pdf)) {
        pdf(out.pdf, 13, 10)
    }
    
    res.df = subset(res.df, !is.na(z) & !is.infinite(z))
    bin.width = median(round(res.df$end - res.df$start + 1))
    ## Pvalue/Qvalue estimation
    if (all(is.na(res.df$z))) 
        return(NULL)
    if (z.th[1] == "sdest") {
        if (min.normal.prop > 0.98) {
            stop("Maximum value accepted for 'min.normal.prop' is 0.98.")
        }
        fdr = fdrtool.quantile(res.df$z, quant.int = seq(min.normal.prop, 0.99, 0.01), 
            plot = !is.null(out.pdf))
        res.df$pv = fdr$pval
        res.df$qv = fdr$qval
        
        ## Remove large aberrations
        aber.large = mergeConsBin.reduce(res.df[which(res.df$qv < 0.05), ], stitch.dist = 10 * 
            bin.width)
        aber.large = subset(aber.large, end - start > 1e+07)
        if (nrow(aber.large) > 0 | !is.null(ref.dist.weight)) {
            if (nrow(aber.large) > 0) {
                aber.gr = with(aber.large, GenomicRanges::GRanges(chr, IRanges::IRanges(start, 
                  end)))
                res.gr = with(res.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, 
                  end)))
                res.aber.large = res.df[GenomicRanges::overlapsAny(res.gr, aber.gr), 
                  ]
                res.df = res.df[!GenomicRanges::overlapsAny(res.gr, aber.gr), ]
            }
            fdr = fdrtool.quantile(res.df$z, quant.int = seq(min.normal.prop, 0.99, 
                0.01), ref.dist.weight = ref.dist.weight, plot = FALSE)
            res.df$pv = fdr$pval
            res.df$qv = fdr$qval
            if (nrow(aber.large) > 0) {
                res.df = rbind(res.df, res.aber.large)
            }
        }
        res.df = res.df[which(res.df$qv < FDR.th), ]
    } else if (z.th[1] == "consbins") {
        res.df = z.thres.cons.bins(res.df, plot = !is.null(out.pdf), pvalues = TRUE)$z.df
    } else if (z.th[1] == "sdest2N") {
        if (min.normal.prop > 0.98) {
            stop("Maximum value accepted for 'min.normal.prop' is 0.98.")
        }
        fdr = fdrtool.quantile.2N(res.df$z, plot = !is.null(out.pdf), min.prop.null = min.normal.prop)
        res.df$pv = fdr$pval
        res.df$qv = fdr$qval
        res.df = res.df[which(res.df$qv < FDR.th), ]
    } else {
        stop("'z.th=': available thresholding approaches are : 'stitch', 'zscores'.")
    }
    
    if (merge.cons.bins[1] != "no") {
        
        if (merge.cons.bins[1] == "stitch") {
            res.df = mergeConsBin.reduce(res.df, stitch.dist = bin.width + 1)
        } else if (merge.cons.bins[1] == "zscores") {
            res.df = mergeConsBin.z(res.df, fdr.th = FDR.th, sd.null = fdr$sigma.est)
        } else {
            stop("'merge.cons.bins=' : available bin merging approaches are : 'stitch', 'zscores'.")
        }
        
        if (nrow(res.df) > 0 & !is.null(out.pdf)) {
            nb.bin.cons = NULL  ## Uglily appease R checks
            print(ggplot2::ggplot(res.df, ggplot2::aes(x = factor(nb.bin.cons))) + 
                ggplot2::geom_histogram() + ggplot2::ylab("number of bins") + ggplot2::xlab("number of consecutive abnormal bins") + 
                ggplot2::theme_bw())
            
            if (any(colnames(res.df) == "fc") & sum(res.df$nb.bin.cons > 2 & res.df$fc < 
                2.5) > 3) {
                print(ggplot2::ggplot(subset(res.df, nb.bin.cons > 2), ggplot2::aes(x = 2 * 
                  fc)) + ggplot2::geom_histogram() + ggplot2::theme_bw() + ggplot2::ylab("number of bins") + 
                  ggplot2::xlab("copy number estimate") + ggplot2::ggtitle("At least 3 consecutive abnormal bins") + 
                  ggplot2::xlim(0, 5))
            }
        }
    }
    
    if (!is.null(out.pdf)) {
        dev.off()
    }
    
    res.df$cn2.dev = abs(res.df$fc - 1)
    
    if (nrow(res.df) > 0 & merge.cons.bins[1] != "no") {
        return(data.frame(sample = samp, res.df, stringsAsFactors = FALSE))
    } else if (any(res.df$qv <= FDR.th, na.rm = TRUE)) {
        return(data.frame(sample = samp, res.df[which(res.df$qv <= FDR.th), ], stringsAsFactors = FALSE))
    } else {
        return(NULL)
    }
} 
