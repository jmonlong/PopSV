##' Quantile normalization of the bin counts. After such normalization, each sample has exactly the same global bin count distribution.
##' @title Quantile normalization
##' @param bc.df a data.frame with 'chr', 'start', 'end' columns and then one column per sample with its bin counts.
##' @param nb.cores the number of cores to use. Default is 1.
##' @return a list with
##' \item{norm.stats}{a data.frame witht some metrics about the normalization of each
##' bin (row) : coverage average and standard deviation; number of outlier reference samples; principal components}
##' \item{bc.norm}{a data.frame, similar to the input 'bc.df', with the normalized bin counts.}
##' @author Jean Monlong
##' @export
quantile.norm <- function(bc.df, nb.cores=1){
  all.samples = setdiff(colnames(bc),c("chr","start","end"))
  rownames(bc) = bins = paste(bc$chr, as.integer(bc$start), as.integer(bc$end), sep="-")
  
  norm.stats = createEmptyDF(c("character", rep("integer",2),rep("numeric",3)), length(bins))
  colnames(norm.stats) = c("chr", "start","end","m","sd","nb.remove")
  bc.norm = createEmptyDF(c("character", rep("integer",2),rep("numeric",length(all.samples))), length(bins))
  colnames(bc.norm) = c("chr", "start","end",all.samples)
  norm.stats$chr = bc.norm$chr = bc$chr
  norm.stats$start = bc.norm$start = bc$start
  norm.stats$end = bc.norm$end = bc$end
  
  bc = as.matrix(bc[,all.samples])
  bc.cov = as.numeric(parallel::mclapply(1:ncol(bc), function(cc) median(bc[,cc], na.rm=TRUE), mc.cores=nb.cores))
  bc = (bc*median(bc.cov)) %*% diag(1/bc.cov)
  bc.sorted = matrix(as.numeric(unlist(parallel::mclapply(1:ncol(bc),function(cc)sort(bc[,cc]),mc.cores=nb.cores))),nrow(bc))
  med = as.numeric(mclapply(1:nrow(bc.sorted), function(rr)median(bc.sorted[rr,], na.rm=TRUE), mc.cores=nb.cores))
  bc.norm[,-(1:3)] = matrix(as.numeric(unlist(paralle::mclapply(1:ncol(bc), function(cc)med[rank(bc[,cc])]),mc.cores=nb.cores)),nrow(bc))
  
  norm.stats[,4:6] = matrix(as.numeric(unlist(parallel::mclapply(1:nrow(bc.norm), function(rr){
    msd = mean.sd.outlierR(as.numeric(bc.norm[rr,all.samples]),1e-6)
    return(c(msd$m,msd$sd,msd$nb.remove))
  },mc.cores=nb.cores))), nrow(bc.norm))
  
  return(list(norm.stats=norm.stats, bc.norm=bc.norm))
}
