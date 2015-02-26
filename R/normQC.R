##' Quality control statistics of bin counts. First, the counts are tested to follow a
##' normal (or Poisson) distribution across samples, in each bin. Then, the randomness
##' of sample ranks are tested. Finally, Z-scores are computed for a subset of the bins
##' and their normality is tested. The second and third test are the most important.
##' Indeed, consistent rankings supports sample-specific technical bias, hence reduced
##' power to detect 'true' abnormal read counts. Non-normal Z-scores will lead to
##' inappropriate fit for the null distribution.
##'
##' Shapiro test is used to test normality of the bin counts across samples. The proportion
##' of bins with non-normal distribution is derived from the Pi0 estimate estimated by
##' package \code{qvalue}. Pi0 is the proportion of pvalues following the null distribution.
##'
##' Goodness of fir from package \code{vcd} is used to test if bin counts follow a Poisson
##' distribution. Again Pi0 estimate from \code{qvalue} package is used to compute the proportion
##' of bin that don't follow Poisson distribution.
##'
##' The randomness of the sample ranks in genomic windows is computed by comparing the position of each sample to the median. If the ranks are random this position should follow a Binomial distribution. For each window we report the number of samples that fail this assumption (Bonferroni corrected P-value<.05). The outputed estimate is the average across all analyzed windows, i.e. the average number of samples with non-random ranks.
##'
##' Z-scores normality is computed by comparing their density distribution and a fitted
##' normal distribution. The estimate represents the proportion of the area under the curve
##' that is unique to the Z-score curve.
##' @title Normalized bin count QC metrics
##' @return a list with
##' \item{prop.non.normal.bin}{proportion of bins with non-normal distribution across samples.}
##' \item{prop.non.poisson.bin}{proportion of bins with non-Poisson distribution across samples.}
##' \item{prop.nonRand.rank}{proportion of bins with non-random ranks.}
##' \item{prop.non.norm.z.mean}{average (across samples) proportion of bins with non-random Z-scores.}
##' \item{prop.non.norm.z.max}{maximum (i.e for worst sample) proportion of bins with non-random Z-scores.}
##' \item{prop.non.norm.z}{proportion of bins with non-random Z-scores, for each sample.}
##' \item{n.subset}{number of bins used for the analysis.}
##' @param bc.df matrix with bin counts (bins x samples). 
##' @param n.subset number of bins to use for the analysis. Default is 10 000. Bins are selected randomly.
##' @param nb.cores the number of cores to use. Default is 1.
##' @author Jean Monlong
##' @export
##' @import magrittr
normQC <- function(bc.df, n.subset = 10000, nb.cores = 1) {
    ## Order by genomic location.
    bc.df = dplyr::arrange(bc.df, chr, start)
    samples = setdiff(colnames(bc.df), c("chr", "start", "end"))
    
    n.subset = min(nrow(bc.df), n.subset)
    
    chrs.t = table(sample(bc.df$chr, n.subset))
    
    test.chr.f <- function(df, win.size = 100) {
        sub.ii = sample.int(nrow(df) - win.size + 1, chrs.t[df$chr[1]])
        res.df = df[sub.ii, c("chr", "start", "end")]
        df = as.matrix(df[, samples])
        ## Bin count normality
        res.df$pv.normal = as.numeric(parallel::mclapply(sub.ii, function(bc.i) {
            bc.i = df[bc.i, ]
            if (length(unique(bc.i)) == 1) return(1)
            return(shapiro.test(bc.i)$p.value)
        }, mc.cores = nb.cores))
        ## Bin count Poisson
        res.df$pv.poisson = as.numeric(parallel::mclapply(sub.ii, function(bc.i) {
            bc.i = df[bc.i, ]
            dump = capture.output({
                res = as.numeric(summary(vcd::goodfit(bc.i, type = "poisson"))[1, 
                  3])
            })
            return(res)
        }, mc.cores = nb.cores))
        ## Ranks randomness
        res.df$nb.rank = as.numeric(parallel::mclapply(sub.ii, function(ii) {
            med.r = apply(df[ii:(ii + win.size - 1), ], 1, median, na.rm = TRUE)
            med.bin = colSums(df[ii:(ii + win.size - 1), ] > med.r)
            pv.bin = pbinom(ifelse(med.bin > 0.5 * length(med.r), length(med.r) - 
                med.bin, med.bin), size = length(med.r), 0.5)
            return(sum(p.adjust(pv.bin, method = "bonf") < 0.05))
        }, mc.cores = nb.cores))
        res.df
    }
    
    chr = . = NULL  ## Uglily appease R checks
    res.df = bc.df %>% dplyr::group_by(chr) %>% dplyr::do(test.chr.f(.))
    
    sub.ii = sample.int(nrow(bc.df), n.subset)
    bc.mat = as.matrix(bc.df[sub.ii, samples])
    ## Z-score distribution
    n.dens = 1000
    z = matrix(as.numeric(unlist(parallel::mclapply(1:nrow(bc.mat), function(bc.i) {
        bc.i = bc.mat[bc.i, ]
        if (any(is.na(bc.i))) {
            return(rep(NA, length(bc.i)))
        }
        bc.i <- as.numeric(bc.i)
        (bc.i - mean(bc.i, na.rm = TRUE))/sd(bc.i, na.rm = TRUE)
    }, mc.cores = nb.cores))), ncol(bc.mat))
    rownames(z) = colnames(bc.mat)
    non.norm.z = as.numeric(parallel::mclapply(1:nrow(z), function(z.samp) {
        z.samp = as.numeric(na.omit(z[z.samp, ]))
        z.norm.est = MASS::fitdistr(z.samp, "normal")$estimate
        f = density(z.samp, n = n.dens, from = min(z.samp), to = max(z.samp))
        fn = dnorm(seq(min(z.samp), max(z.samp), length.out = n.dens), z.norm.est[1], 
            z.norm.est[2])
        dis.prop = sum(abs(f$y - fn))/(2 * sum(fn))
        return(dis.prop)
    }, mc.cores = nb.cores))
    
    ## PCA dispersion
    pca.o = prcomp(t(bc.mat))
    pca.d = as.matrix(dist(pca.o$x[, 1:2]))
    pca.dmm = mean(pca.d)/median(pca.d)
    
    qv.normal = qvalue::qvalue(res.df$pv.normal)
    if (mean(res.df$pv.poisson < 0.05, na.rm=TRUE) > 0.9) {
        qv.poisson = list(pi0 = 0)
    } else {
        qv.poissson = qvalue::qvalue(res.df$pv.poisson)
    }
    ## 
    return(list(prop.non.normal.bin = 1 - qv.normal$pi0, prop.non.poisson.bin = 1 - 
        qv.poisson$pi0, nb.nonRand.rank = mean(res.df$nb.rank), prop.nonNorm.z.mean = mean(non.norm.z), 
        prop.nonNorm.z.max = max(non.norm.z), prop.nonNorm.z = non.norm.z, pca.dmm = pca.dmm, 
        n.subset = n.subset))
} 
