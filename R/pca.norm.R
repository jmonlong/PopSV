##' Bin counts are normalized by regressing out the effect of the first Principal
##' Components. Beforehands, the average coverage is normalized. Then PC are computed
##' using \code{prcomp} function and regressed out using linear regression.
##' @title PCA-based normalization of bin counts
##' @param bc bin count matrix or data.frame (bins x samples)
##' @param nb.pcs the number of Principal Components to include in the regression model.
##' @param ref.samples the names of the reference samples on which the PCA is performed.
##' By default all samples are used.
##' @return a matrix with the normalized bin counts (bin x sample)
##' @author Jean Monlong
##' @export
pca.norm <- function(bc, nb.pcs=3, ref.samples=NULL){
  bc.cov = colMeans(bc, na.rm=TRUE)
  bc = t(t(bc)*mean(bc.cov)/bc.cov)
  if(is.null(ref.samples)) ref.samples = colnames(bc)
  bc = bc[!is.na(bc[,1]),]
  pca.o = prcomp(bc[, ref.samples])
  pca.o = pca.o$x[,1:nb.pcs]
  colnames(pca.o) = paste("pc",1:nb.pcs,sep="")
  reg.form = paste("bc.s ~ pc", paste(1:nb.pcs, collapse=" + pc"),sep="")
  apply(bc, 2, function(bc.s){
    lm.o = lm(reg.form, data=data.frame(bc.s=bc.s, pca.o))
    return(lm.o$coefficients[1]+lm.o$residuals)
  })
}
