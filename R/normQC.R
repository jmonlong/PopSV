##' Quality control statistics of bin counts. First, the counts are tested to follow a
##' normal (or Poisson) distribution across samples, in each bin. Then, the randomness
##' of sample ranks are tested. Finally, Z-scores are computed for a subset of the bins
##' and their normality is tested. The second and third test are the most important.
##' Indeed, consistent rankings supports sample-specific technical bias, hence reduced
##' power to detect "true" abnormal read counts. Non-normal Z-scores will lead to
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
##' The randomness of the sample ranks are tested using a Chi-Square test. To have an
##' idea of the proportion of the bins with abnormal ranking, bins are randomly divided
##' into groups and the test performed in each group. The proportion of bin with non-random
##' ranks is approximated by the proportion of groups which failed the Chi-Square test.
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
##' @author Jean Monlong
##' @export
normQC <- function(bc.df, n.subset=1e4){
    ## Order by genomic location.
    bc.df = arrange(bc.df, chr, start)
    samples = setdiff(colnames(bc.df), c("chr","start","end"))

    n.subset = min(nrow(bc.mat), n.subset)

    chrs.t = table(sample(bc.df$chr,n.subset))

    test.chr.f <- function(df, win.size=100){
        sub.ii = sample.int(nrow(df)-win.size+1, chrs.t[df$chr[1]])
        res.df = df[, c("chr","start","end")]
        df = as.matrix(df[,samples])
        ## Bin count normality
        res.df$pv.normal = as.numeric(apply(df[sub.ii,], 1, function(bc.i){
            if(length(unique(bc.i))==1) return(1)
            return(shapiro.test(bc.i)$p.value)
        }))
        
        ## Bin count Poisson
        res.df$pv.poisson = as.numeric(apply(df[sub.ii,], 1, function(bc.i){
            return(as.numeric(summary(vcd::goodfit(bc.i, type="poisson"))[1,3]))
        }))

        ## Ranks randomness
        res.df$pv.rank = sapply(sub.ii, function(ii){
            rank.mat = apply(df[ii:(ii+win.size-1),], 1, function(bc.i){
                if(any(!is.na(bc.i))){
                    rank(bc.i,ties.method="random")
                } else {
                    rep(NA,length(bc.i))
                }
            })
            rk.t = apply(rank.mat, 1, function(rks){
                table(factor(rks,levels=1:nrow(rank.mat)))
            })
            return(suppressWarnings(chisq.test(rk.t)$p.value))
        })

        res.df
    }

    res.df = bc.df %>% group_by(chr) %>% do(test.chr.f)
    
    sub.ii = sample.int(nrow(bc.df), n.subset)
    bc.mat = as.matrix(bc.df[,samples])
    ## Z-score distribution
    n.dens = 1000
    z = apply(bc.mat,1,function(bc.i){
        if(any(is.na(bc.i))){
            return(rep(NA, length(bc.i)))
        }
        bc.i <- as.numeric(bc.i)
        (bc.i-mean(bc.i,na.rm=TRUE))/sd(bc.i,na.rm=TRUE)
    })
    rownames(z) = colnames(bc.mat)
    non.norm.z = apply(z, 1, function(z.samp){
        z.samp = as.numeric(na.omit(z.samp))
        z.norm.est = MASS::fitdistr(z.samp,"normal")$estimate
        f = density(z.samp,n=n.dens,from=min(z.samp),to=max(z.samp))
        fn = dnorm(seq(min(z.samp),max(z.samp),length.out=n.dens),z.norm.est[1],z.norm.est[2])
        dis.prop = sum(abs(f$y-fn)) / (2 * sum(fn))
        return(dis.prop)
    })

    ## PCA entropy
    pca.o = prcomp(t(bc.mat))
    pca.d = as.matrix(dist(pca.o$x[,1:2]))
    pca.dmm = mean(pca.d) / median(pca.d)
    bc.shuf = apply(bc.mat, 1, sample)
    pca.s = prcomp(bc.shuf)
    pca.d.s = as.matrix(dist(pca.s$x[,1:2]))
    pca.dmm.s = mean(pca.d.s) / median(pca.d.s)
    
    qv.normal = qvalue::qvalue(res.df$pv.normal)    
    qv.poissson = qvalue::qvalue(res.df$pv.poisson)    
    ## 
    return(list(prop.non.normal.bin = 1-qv.normal$pi0,
                prop.non.poisson.bin = 1-qv.poisson$pi0,
                prop.nonRand.rank = mean(res.df$pv.rank<=.05),
                prop.nonNorm.z.mean = mean(non.norm.z),
                prop.nonNorm.z.max = max(non.norm.z),
                prop.nonNorm.z = non.norm.z,
                pca.dmm = pca.dmm,
                pca.dmm.s = pca.dmm.s,
                n.subset = n.subset))
}
