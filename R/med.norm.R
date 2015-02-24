##' Global median normalization of bin counts. After such normalization, each sample has the same median bin count. 
##' @title Global median normalization
##' @param bc.df a data.frame with 'chr', 'start', 'end' columns and then one column per sample with its bin counts.
##' @param nb.cores the number of cores to use. Default is 1.
##' @param norm.stats.comp Should some statistics on the normalized bin count be computed (mean, sd, outliers). Default is TRUE.
##' @return a list with
##' \item{norm.stats}{a data.frame witht some metrics about the normalization of each
##' bin (row) : coverage average and standard deviation; number of outlier samples}
##' \item{bc.norm}{a data.frame, similar to the input 'bc.df', with the normalized bin counts.}
##' @author Jean Monlong
##' @export
med.norm <- function(bc.df, nb.cores=1, norm.stats.comp=TRUE){
  all.samples = setdiff(colnames(bc.df),c("chr","start","end"))
  rownames(bc.df) = bins = paste(bc.df$chr, as.integer(bc.df$start), as.integer(bc.df$end), sep="-")
  
  norm.stats = createEmptyDF(c("character", rep("integer",2),rep("numeric",3)), length(bins))
  colnames(norm.stats) = c("chr", "start","end","m","sd","nb.remove")
  bc.norm = createEmptyDF(c("character", rep("integer",2),rep("numeric",length(all.samples))), length(bins))
  colnames(bc.norm) = c("chr", "start","end",all.samples)
  norm.stats$chr = bc.norm$chr = bc.df$chr
  norm.stats$start = bc.norm$start = bc.df$start
  norm.stats$end = bc.norm$end = bc.df$end
  
  bc.mat = as.matrix(bc.df[,all.samples])
  bc.cov = as.numeric(parallel::mclapply(1:ncol(bc.mat), function(cc) median(bc.mat[,cc], na.rm=TRUE), mc.cores=nb.cores))
  bc.norm[,-(1:3)] = (bc.mat*median(bc.cov)) %*% diag(1/bc.cov)

  if(norm.stats.comp){
    norm.stats[,4:6] = matrix(as.numeric(unlist(parallel::mclapply(1:nrow(bc.norm), function(rr){
      msd = mean.sd.outlierR(as.numeric(bc.norm[rr,all.samples]),1e-6)
      return(c(msd$m,msd$sd,msd$nb.remove))
    },mc.cores=nb.cores))), nrow(bc.norm))
  } else {
    norm.stats.comp = NULL
  }
  
  return(list(norm.stats=norm.stats, bc.norm=bc.norm))
}
