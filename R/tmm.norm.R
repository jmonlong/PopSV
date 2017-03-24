##' TMM normalization of bin counts. After such normalization, the majority of bins have the same trimmed average across samples.
##' @title Trimmed-Mean M normalization
##' @param bc.df a data.frame with 'chr', 'start', 'end' columns and then one column per sample with its bin counts.
##' @param cont.sample the sample to use as baseline for the pairwise normalization.
##' All the samples will be normalized to it.
##' @param trim.level How much of the M values should be trimmed (on each side of the distribution). Default is 0.3 (30\%).
##' @param nb.cores the number of cores to use. Default is 1.
##' @param norm.stats.comp Should some statistics on the normalized bin count be computed (mean, sd, outliers). Default is TRUE.
##' @return a list with
##' \item{norm.stats}{a data.frame witht some metrics about the normalization of each
##' bin (row) : coverage average and standard deviation; number of outlier samples}
##' \item{bc.norm}{a data.frame, similar to the input 'bc.df', with the normalized bin counts.}
##' @author Jean Monlong
##' @export
tmm.norm <- function(bc.df, cont.sample, trim.level=.3, nb.cores = 1, norm.stats.comp = TRUE) {
    all.samples = setdiff(colnames(bc.df), c("chr", "start", "end", "chr2", "start2"))

    if(all(all.samples != cont.sample)){
        stop("'cont.sample' sample is not in 'bc.df'")
    }

    bc.norm = createEmptyDF(rep("numeric", length(all.samples)),nrow(bc.df))
    colnames(bc.norm) = all.samples
    cols = intersect(colnames(bc.df), c("chr", "start", "end", "chr2", "start2"))
    bc.norm = cbind(bc.df[,cols], bc.norm)

    bc.mat = as.matrix(bc.df[, all.samples])
    bc.cont = bc.mat[,cont.sample]
    norm.factor = as.numeric(parallel::mclapply(1:ncol(bc.mat), function(cc) {
                                                    bc.samp = bc.mat[,cc]
                                                    no00 = which(bc.samp!=0 & bc.cont!=0)
                                                    exp(mean(log(bc.cont[no00]/bc.mat[no00,cc]), trim=trim.level))
                                                }, mc.cores = nb.cores))

    bc.norm[, all.samples] = round(bc.mat %*% diag(norm.factor), 2)

    if (norm.stats.comp) {
      norm.stats = matrix(as.numeric(unlist(parallel::mclapply(1:nrow(bc.norm),
          function(rr) {
              msd = mean.sd.outlierR(as.numeric(bc.norm[rr, all.samples]))
              return(c(msd$m, msd$sd, msd$nb.remove))
          }, mc.cores = nb.cores))), nrow(bc.norm), byrow=TRUE)
      norm.stats = as.data.frame(norm.stats)
      colnames(norm.stats) = c("m", "sd", "nb.remove")
      norm.stats = cbind(bc.df[,cols], norm.stats)
    } else {
        norm.stats = NULL
    }

    return(list(norm.stats = norm.stats, bc.norm = bc.norm))
}
