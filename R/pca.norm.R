##' Bin counts are normalized by regressing out the effect of the first Principal
##' Components. Beforehands, the average coverage is normalized. Then PC are computed
##' using \code{prcomp} function and regressed out using linear regression.
##' @title PCA-based normalization of bin counts
##' @param bc bin count matrix or data.frame (bins x samples)
##' @param nb.pcs the number of Principal Components to include in the regression model.
##' @param ref.samples the names of the reference samples on which the PCA is performed.
##' By default all samples are used.
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @return a list with
##' \item{norm.stats}{a data.frame witht some metrics about the normalization of each
##' bin (row) : coverage average and standard deviation; number of outlier reference samples; principal components}
##' \item{bc.norm}{a matrix with the normalized bin counts (bin x sample).}
##' \item{z}{a matrix with the Z-scores for each bin and sample (bin x sample).}
##' \item{fc}{a matrix with the copy-number coefficients for each bin and sample (
##' bin x sample).}
##' \item{z.poisson}{a backup of the input parameters.}
##' @author Jean Monlong
##' @export
pca.norm <- function(bc, nb.pcs=3, ref.samples=NULL, z.poisson=FALSE){
    all.samples = setdiff(colnames(bc),c("chr","start","end"))
    if(is.null(ref.samples)) ref.samples = all.samples
    rownames(bc) = bins = paste(bc$chr, as.integer(bc$start), sep="-")
    norm.stats = createEmptyDF(c("character", rep("integer",2),rep("numeric",3+nb.pcs)), length(bins))
    colnames(norm.stats) = c("chr", "start","end","m","sd","nb.remove", paste0("PC",1:nb.pcs))
    bc.norm = createEmptyDF(c("character", rep("integer",2),rep("numeric",length(all.samples))), length(bins))
    z = createEmptyDF(c("character", rep("integer",2),rep("numeric",length(all.samples))), length(bins))
    fc = createEmptyDF(c("character", rep("integer",2),rep("numeric",length(all.samples))), length(bins))
    colnames(bc.norm) = colnames(z) = colnames(fc) = c("chr", "start","end",all.samples)
    norm.stats$chr = bc.norm$chr = z$chr = fc$chr = bc$chr
    norm.stats$start = bc.norm$start = z$start = fc$start = bc$start
    norm.stats$end = bc.norm$end = z$end = fc$end = bc$end
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
    bc = bc[,all.samples]
    bc.cov = colMeans(bc, na.rm=TRUE)
    bc = t(t(bc)*mean(bc.cov)/bc.cov)
    if(any(is.na(bc))){
        bc[is.na(bc)] = 0
    }
    pca.o = prcomp(bc[,ref.samples])
    pca.o = pca.o$x[,1:nb.pcs]
    colnames(pca.o) = paste("pc",1:nb.pcs,sep="")
    reg.form = paste("bc.s ~ pc", paste(1:nb.pcs, collapse=" + pc"),sep="")
    bc.norm[, all.samples] = apply(bc, 2, function(bc.s){
        lm.o = lm(reg.form, data=data.frame(bc.s=bc.s, pca.o))
        return(bc.s*lm.o$coefficients[1]/predict(lm.o))
    })
    for(bin.ii in 1:length(bins)){
        msd = mean.sd.outlierR(as.numeric(bc.norm[bin.ii,ref.samples]),1e-6)
        norm.stats[bin.ii,4:6]=c(msd$m,msd$sd,msd$nb.remove)
        z[bin.ii,-(1:3)] = z.comp(bc.norm[bin.ii,],msd$m,msd$sd)
        fc[bin.ii,-(1:3)] = bc.norm[bin.ii,]/msd$m
    }
    norm.stats[,-(1:6)] = pca.o
    list(norm.stats=norm.stats, bc.norm=bc.norm,z=z, fc=fc, z.poisson=z.poisson)
}
