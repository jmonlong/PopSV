##' Bin counts are normalized by regressing out the effect of the first Principal
##' Components. Beforehands, the average coverage is normalized. Then PC are computed
##' using \code{prcomp} function and regressed out using linear regression.
##' @title PCA-based normalization of bin counts
##' @param bc bin count matrix or data.frame (bins x samples)
##' @param nb.pcs the number of Principal Components to include in the regression model.
##' @return a list with
##' \item{norm.stats}{a data.frame witht some metrics about the normalization of each
##' bin (row) : coverage average and standard deviation; number of outlier reference samples; principal components}
##' \item{bc.norm}{a matrix with the normalized bin counts (bin x sample).}
##' @author Jean Monlong
##' @export
pca.norm <- function(bc, nb.pcs=3, nb.cores=1){
  
  all.samples = setdiff(colnames(bc),c("chr","start","end"))
  rownames(bc) = bins = paste(bc$chr, as.integer(bc$start), as.integer(bc$end), sep="-")
  
  norm.stats = createEmptyDF(c("character", rep("integer",2),rep("numeric",3+nb.pcs)), length(bins))
  colnames(norm.stats) = c("chr", "start","end","m","sd","nb.remove", paste0("PC",1:nb.pcs))
  bc.norm = createEmptyDF(c("character", rep("integer",2),rep("numeric",length(all.samples))), length(bins))
  colnames(bc.norm) = c("chr", "start","end",all.samples)
  norm.stats$chr = bc.norm$chr = bc$chr
  norm.stats$start = bc.norm$start = bc$start
  norm.stats$end = bc.norm$end = bc$end
  
  bc = as.matrix(bc[,all.samples])
  bc.cov = as.numeric(parallel::mclapply(1:ncol(bc), function(cc) median(bc[,cc], na.rm=TRUE), mc.cores=nb.cores))
  bc = (bc*median(bc.cov)) %*% diag(1/bc.cov)
  if(any(is.na(bc))){
    bc[is.na(bc)] = 0
  }

  pca.o = prcomp(bc)
  pca.o = pca.o$x[,1:nb.pcs]
  colnames(pca.o) = paste("pc",1:nb.pcs,sep="")
  reg.form = paste("bc.s ~ pc", paste(1:nb.pcs, collapse=" + pc"),sep="")
  bc.norm[, all.samples] = matrix(as.numeric(unlist(mclapply(1:ncol(bc), function(cc){
    lm.o = glm(reg.form, data=data.frame(bc.s=bc[,cc], pca.o)) ##, family=poisson)
    return(bc.s*lm.o$coefficients[1]/predict(lm.o))
  },mc.cores=nb.cores))), nrow(bc))

  norm.stats[,4:6] = matrix(as.numeric(unlist(mclapply(1:nrow(bc.norm), function(rr){
    msd = mean.sd.outlierR(as.numeric(bc.norm[rr,all.samples]),1e-6)
    return(c(msd$m,msd$sd,msd$nb.remove))
  },mc.cores=nb.cores))), nrow(bc.norm))
  norm.stats[,-(1:6)] = pca.o

  return(list(norm.stats=norm.stats, bc.norm=bc.norm))
}
