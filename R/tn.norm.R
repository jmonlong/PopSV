##' Bin counts are normalized one bin at a time, using a subset of the bins that look
##' similar across the reference samples.
##'
##' The Z-score is computed by substracting the bin count by the average bin count
##' across the reference samples and dividing by their standard deviation. If
##' 'z.poisson' is TRUE, a score using Poisson distribution is also computed, using
##' the average bin count as an estimator of the lambda. Then the score with the lowest
##' absolute value is kept. This hybrid Z-score is to be used when some regions have low
##' coverage where it is more robust to use Poisson assumptions.
##' @title Targeted-normalization of bin counts
##' @param bc a matrix or data.frame with the bin counts (bin x sample).
##' @param cont.sample the sample to use as baseline for the pairwise normalization.
##' All the samples will be normalized to it.
##' @param ref.samples a vector with the names of the samples to be used as reference.
##' @param nb.support.bins the number of bins to use for the normalization.
##' @param bins a vector the names of the bins to normalize. If NULL (default), all
##' bins are normalized.
##' @param save.support.bins if TRUE (default) the bins used for the normalization are
##' saved in the output object 'norm.stats'.
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @return a list with
##' \item{norm.stats}{a data.frame witht some metrics about the normalization of each
##' bin (row) : correlation with worst supporting bin; coverage average and standard
##' deviation; number of outlier reference samples; supporting bins.}
##' \item{bc.norm}{a matrix with the normalized bin counts (bin x sample).}
##' \item{z}{a matrix with the Z-scores for each bin and sample (bin x sample).}
##' \item{cn.coeff}{a matrix with the copy-number coefficients for each bin and sample (
##' bin x sample).}
##' \item{nb.support.bins, cont.sample, z.poisson}{a backup of the input parameters.}
##' @author Jean Monlong
##' @export
tn.norm <- function(bc,cont.sample,ref.samples,nb.support.bins=1e3,bins=NULL,save.support.bins=TRUE, z.poisson=FALSE){
    all.samples = colnames(bc)
    ref.samples.ii = which(colnames(bc) %in% ref.samples)
    bc = t(bc)
    if(is.null(bins)) bins = colnames(bc)
    if(z.poisson){
        z.comp <- function(x, mean.c, sd.c){
            z.n = (x-mean.c)/sd.c
            z.p = qnorm(ppois(x, mean.c))
            n.ii = abs(z.n) < abs(z.p)
            z.p[n.ii] = z.n[n.ii]
            z.p
        }
    } else {
        z.comp <- function(x, mean.c, sd.c){(x-mean.c)/sd.c}
    }
    if(save.support.bins) {
        norm.stats = createEmptyDF(c("character",rep("numeric",4),rep("character",nb.support.bins)), length(bins))
        colnames(norm.stats) = c("bin","d.max","mean","sd","nb.remove",paste("b",1:nb.support.bins,sep=""))
    } else {
        norm.stats = createEmptyDF(c("character",rep("numeric",4)), length(bins))
        colnames(norm.stats) = c("bin","d.max","mean","sd","nb.remove")
    }
    bc.norm = createEmptyDF(c("character",rep("numeric",nrow(bc))), length(bins))
    z = createEmptyDF(c("character",rep("numeric",nrow(bc))), length(bins))
    cn.coeff = createEmptyDF(c("character",rep("numeric",nrow(bc))), length(bins))
    colnames(bc.norm) = colnames(z) = colnames(cn.coeff) = c("bin",all.samples)
    norm.stats$bin = bc.norm$bin = z$bin = cn.coeff$bin = bins
    for(bin.ii in 1:length(bins)){
        bin = bins[bin.ii]
        bc.i = bc[,bin]
        if(any(!is.na(bc.i) & bc.i!=0)) {
            if(sum(as.numeric(bc.i[ref.samples.ii])>0)<5){
                d.o.i = c(which(colnames(bc)==bin),sample(1:ncol(bc), nb.support.bins-1))
                d.max=-1
            } else {
                d.i = 1-as.numeric(cor(as.numeric(bc.i[ref.samples.ii]),bc[ref.samples.ii,],use="pairwise.complete.obs"))
                d.o.i = order(d.i)[1:nb.support.bins]
                d.max = as.numeric(d.i[d.o.i[nb.support.bins]])
            }
            if(!is.infinite(d.max) & !is.na(d.max)) {
                bc.g = t(bc[,d.o.i])
                bin.for.norm = colnames(bc)[d.o.i]
                norm.coeff = rep(NA,ncol(bc.g))
                norm.coeff[ref.samples.ii] = norm.tm.opt(bc.g[,ref.samples.ii],ref.col=bc.g[,cont.sample])
                msd = mean.sd.outlierR(bc.g[1,ref.samples.ii] * norm.coeff[ref.samples.ii],1e-6)
                if(length(norm.coeff) > length(ref.samples.ii)){
                    norm.coeff[-ref.samples.ii] = norm.tm.opt(bc.g[,-ref.samples.ii],ref.col=bc.g[,cont.sample],bc.mean.norm=msd$mean)
                }
                bc.t = bc.g[1,] * norm.coeff
                if(any(!is.na(bc.t))){
                    norm.stats[bin.ii,2:5]=c(d.max,msd$mean,msd$sd,msd$nb.remove)
                    if(save.support.bins) norm.stats[bin.ii,6:ncol(norm.stats)] = bin.for.norm
                    bc.norm[bin.ii,-1] = bc.t
                    z[bin.ii,-1] = z.comp(bc.t,msd$mean,msd$sd)
                    cn.coeff[bin.ii,-1] = bc.t/msd$mean

                }
            }
        }
    }
    list(norm.stats=norm.stats, bc.norm=bc.norm,z=z, cn.coeff=cn.coeff, nb.support.bins=nb.support.bins, cont.sample=cont.sample, z.poisson=z.poisson)
}
