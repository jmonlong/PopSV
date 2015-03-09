##' Z-score computation from bin count and/or mean/sd metrics on the reference samples
##'
##' @title Z-score computation
##' @param files.df a data.frame with the file paths.
##' @param samples the samples to analyze.
##' @param msd.f the path to the file with mean/sd bin count on reference samples. 
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @param col the column in 'files.df' that define the bin count file path.
##' @param nb.cores the number of cores to use.
##' @return a list with
##' \item{z}{a data.frame with the Z-scores for each bin and sample (bin x sample).}
##' \item{fc}{a data.frame with the fold-change compared to the average bin count in
##' the reference samples for each bin and sample (bin x sample).}
##' \item{msd}{the mean, standard deviation and number of removed outlier samples in each bin.}
##' @author Jean Monlong
##' @export
z.comp <- function(files.df, samples, msd.f = NULL, z.poisson = FALSE, col = "bc.gc.norm.gz", 
    nb.cores = 1) {
    
    if (z.poisson) {
        z.comp.f <- function(x, mean.c, sd.c) {
            z.n = (x - mean.c)/sd.c
            z.p = qnorm(ppois(x, mean.c))
            n.ii = abs(z.n) < abs(z.p)
            z.p[which(n.ii)] = z.n[which(n.ii)]
            z.p
        }
    } else {
        z.comp.f <- function(x, mean.c, sd.c) {
            (x - mean.c)/sd.c
        }
    }
        
    if (is.data.frame(files.df)) {
        bc.1 = data.table::fread(subset(files.df, sample == samples[1])[, col], header = TRUE)
        bc.l = parallel::mclapply(subset(files.df, sample %in% samples)[, col], function(fi) {
          data.table::fread(fi, header = TRUE)[, 4, with = FALSE]
        }, mc.cores = nb.cores)
        bc.l = matrix(unlist(bc.l), ncol=length(bc.l))
    } else {
        bc.l = data.table::fread(files.df, header = TRUE)
        bc.1 = bc.l[, 1:3, with = FALSE]
        bc.l = as.matrix(bc.l[, samples, with = FALSE])
    }
    
    if (is.null(msd.f)) {
        msd = parallel::mclapply(1:nrow(bc.l), function(rr) unlist(mean.sd.outlierR(bc.l[rr,], pv.max.ol = 1e-06)), mc.cores=nb.cores)
        msd = matrix(unlist(msd), nrow=3)
    } else {
        msd = as.data.frame(data.table::fread(msd.f, select = 1:5, header = TRUE))
        msd.col.ids = 1:nrow(msd)
        names(msd.col.ids) = paste(msd$chr, as.integer(msd$start), as.integer(msd$end), sep = "-")
        msd = t(as.matrix(msd[, -(1:3)]))
        msd = msd[, msd.col.ids[paste(bc.1$chr, as.integer(bc.1$start), as.integer(bc.1$end), sep = "-")]]
    }
    
    z = parallel::mclapply(1:ncol(bc.l), function(cc) z.comp.f(bc.l[,cc], mean.c = msd[1, ], sd.c = msd[2, ]), mc.cores=nb.cores)
    z = matrix(unlist(z), ncol=length(z))
    fc = bc.l/msd[1, ]
    colnames(z) = colnames(fc) = samples
    z = data.frame(bc.1[, 1:3, with = FALSE], z)
    fc = data.frame(bc.1[, 1:3, with = FALSE], fc)
    if (is.null(msd.f)) {
        msd = data.frame(bc.1[, 1:3], t(msd))
    } else {
        msd = NULL
    }
    return(list(z = z, fc = fc, msd = msd, z.poisson = z.poisson))
}

## To test: one sample only; chunks; msd input or not 
