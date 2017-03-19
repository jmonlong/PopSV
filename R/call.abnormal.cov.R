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
##' @param sub.z if non-NULL the number of bins in a sub-segment for Z-score null distribution estimation. Default is NULL. If highly rearranged genomes (cancer), try '1e4'.
##' @param outfile.pv if non-NULL, the name of the file to write all the Pvalues (for all bins). Used in some analysis (e.g. annotate.with.parents).
##' @return a data.frame with columns
##' \item{chr, start, end}{the genomic region definition.}
##' \item{z}{the Z-score.}
##' \item{pv, qv}{the P-value and Q-value(~FDR).}
##' \item{fc}{the copy number estimate (if 'fc' was not NULL).}
##' \item{nb.bin.cons}{the number of consecutive bins (if the bins were merged, i.e. ' 'merge.cons.bins!='no'').}
##' \item{cn2.dev}{Copy number deviation from the reference.}
##' @author Jean Monlong
##' @export
call.abnormal.cov <- function(files.df, samp, out.pdf = NULL, FDR.th = 0.05, merge.cons.bins = c("stitch", "zscores", "cbs", "no"), stitch.dist=NULL, z.th = c("sdest", "consbins", "sdest2N"), norm.stats = NULL, min.normal.prop = 0.9, aneu.chrs = NULL, gc.df=NULL, sub.z=NULL, outfile.pv=NULL) {

    if(!is.data.frame(files.df) & is.data.frame(files.df$z)){
        res.df = files.df$z
        if(all(samp != colnames(res.df))){
            stop(samp, " is not in the z-score data.frame.")
        }
        res.df = res.df[,c(colnames(res.df)[1:3], samp)]
        colnames(res.df)[4] = "z"
        if(is.data.frame(files.df$fc) & samp %in% colnames(files.df$fc)){
            res.df$fc = files.df$fc[,samp]
        }
    } else {
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
    }

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
        if(is.null(sub.z)){
            fdr = fdrtool.quantile(res.df$z, quant.int = seq(min.normal.prop, 0.99, 0.01))
        } else {
            ## Get 100 sub-segments of 'sub.z' consecutive bins and use the 20 with median closest to 0 (to focus on "normal" genome and avoid large aberrations)
            res.df = res.df[order(res.df$chr, res.df$start),]
            subz.l = lapply(round(stats::runif(100,1,nrow(res.df)-sub.z-1)), function(ss) res.df[ss:(ss+sub.z),])
            subz.med = unlist(lapply(subz.l, function(df) stats::median(df$z, na.rm=TRUE)))
            subz.z = unlist(lapply(subz.l[order(abs(subz.med))[1:20]], function(df)df$z))
            fdr = fdrtool.quantile(subz.z, quant.int = seq(min.normal.prop, 0.99, 0.01))
        }
    } else if (z.th[1] == "consbins") {
        fdr = z.thres.cons.bins(res.df)
    } else if (z.th[1] == "sdest2N") {
        fdr = fdrtool.quantile.2N(res.df$z)
    } else {
        stop("'z.th=': available thresholding approaches are : 'stitch', 'zscores'.")
    }

    res.df$qv = res.df$pv = NA
    res.df$pv[which(res.df$z > 0)] = 2 * stats::pnorm(-abs(res.df$z[which(res.df$z > 0)]), 0, fdr$sigma.est.dup)
    res.df$pv[which(res.df$z < 0)] = 2 * stats::pnorm(-abs(res.df$z[which(res.df$z < 0)]), 0, fdr$sigma.est.del)
    if (any(res.df$pv == 0, na.rm = TRUE)){
        res.df$pv[which(res.df$pv == 0)] = .Machine$double.xmin
    }
    res.df$qv = stats::p.adjust(res.df$pv, method = "fdr")

    if(!is.null(outfile.pv)){
        pv.df = res.df[,c("chr","start","end","pv","qv")]
        pv.df = pv.df[order(pv.df$chr, pv.df$start),]
        utils::write.table(pv.df, outfile.pv, quote=FALSE, row.names=FALSE, sep="\t")
        comp.index.files(outfile.pv)
        rm(pv.df)
    }

    if (!is.null(out.pdf) & any(!is.na(res.df$pv))) {
        pv = qv = ..density.. = y = z = NULL  ## Uglily appease R checks

        z.lim = c(-fdr$sigma.est.del, fdr$sigma.est.dup)*ifelse(mean(res.df$pv<.01)>.1,8,5)
        null.df = data.frame(y=c(stats::dnorm(seq(z.lim[1],0,.05),0,fdr$sigma.est.del),stats::dnorm(seq(0,z.lim[2],.05),0,fdr$sigma.est.dup)), z=c(seq(z.lim[1],0,.05),seq(0,z.lim[2],.05)))
        null.df$y = null.df$y * mean(res.df$z> -4*fdr$sigma.est.del & res.df$z<4*fdr$sigma.est.dup)

        print(ggplot2::ggplot(res.df, ggplot2::aes(x = z)) +
                  ggplot2::geom_histogram(ggplot2::aes(y=..density..), bins=30, na.rm=TRUE) +
                      ggplot2::xlab("Z-score") + ggplot2::ylab("number of bins") + ggplot2::theme_bw() +
                          ggplot2::geom_line(ggplot2::aes(y=y), data=null.df, linetype=2, colour="red") +
                              ggplot2::xlim(z.lim))

        print(ggplot2::ggplot(res.df, ggplot2::aes(x = pv, fill=cut(qv, breaks = c(-Inf, 0.001, 0.01, 0.5, 0.1,1)))) + ggplot2::geom_histogram(bins=30, na.rm=TRUE) +
                  ggplot2::xlab("P-value") + ggplot2::xlim(-0.2, 1) + ggplot2::ylab("number of bins") +
                      ggplot2::scale_fill_hue(name="Q-value") +
                          ggplot2::theme_bw() + ggplot2::theme(legend.position="bottom"))
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
            res.df = mergeConsBin.cbs(res.df, pv.th = FDR.th, stitch.dist = stitch.dist)
        } else {
            stop("'merge.cons.bins=' : available bin merging approaches are : 'stitch', 'zscores'.")
        }

        if (nrow(res.df) > 0 & !is.null(out.pdf)) {
            nb.bin.cons = NULL  ## Uglily appease R checks
            print(ggplot2::ggplot(res.df, ggplot2::aes(x = factor(nb.bin.cons))) +
                      ggplot2::geom_bar(na.rm=TRUE) + ggplot2::ylab("number of bins") + ggplot2::xlab("number of consecutive abnormal bins") +
                          ggplot2::theme_bw())

            if (any(colnames(res.df) == "fc") & sum(res.df$nb.bin.cons > 2 & res.df$fc < 2.5) > 3) {
                print(ggplot2::ggplot(subset(res.df, nb.bin.cons > 2), ggplot2::aes(x = 2 * fc)) + ggplot2::geom_histogram(bins=30, na.rm=TRUE) + ggplot2::theme_bw() + ggplot2::ylab("number of bins") + ggplot2::xlab("copy number estimate") + ggplot2::ggtitle("At least 3 consecutive abnormal bins") + ggplot2::xlim(0, 5))
            }
        }
    }

    if (!is.null(out.pdf)) {
        grDevices::dev.off()
    }

    ## Deviation from copy number 2
    res.df$cn2.dev = round(2*abs(res.df$fc - 1), 4)
    res.df$cn = round(2*res.df$fc)
    if(nrow(res.df)>0){
        res.df$prop.single.bin = round(mean(res.df$nb.bin.cons==1), 3)
        res.df$pv = signif(res.df$pv, 4)
        res.df$qv = signif(res.df$qv, 4)
        res.df$z = round(res.df$z, 2)
        res.df$fc = round(res.df$fc, 5)
    } else {
        res.df$prop.single.bin = res.df$pv = res.df$qv = res.df$z = res.df$fc = numeric(0)
    }

    if (nrow(res.df) > 0 & merge.cons.bins[1] != "no") {
        return(data.frame(sample = samp, res.df, stringsAsFactors = FALSE))
    } else if (any(res.df$qv <= FDR.th, na.rm = TRUE)) {
        return(data.frame(sample = samp, res.df[which(res.df$qv <= FDR.th), ], stringsAsFactors = FALSE))
    } else {
        return(data.frame())
    }
}
