##' This function retrieves bin count for all samples in the input file information
##' data.frame. It was designed to count reads in a small amount of bins (a subset
##' of what will be used for the full analysis). This counts can be used later to
##' compare a potentially large number samples and guide the selection of samples
##' to consider in the analysis. 
##' @title Counts reads across samples in a small number of bins
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bam' are required.
##' @param bins.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param nb.cores number of cores to use. If higher than 1, \code{parallel}
##' package is used to parallelized the counting.
##' @param col.files if non-NULL, the name of the column in 'files.df' with the filenames where the information should be retrieved from.
##' @param ... extra paramaters for 'bin.bam' function
##' @return a data.frame with the bin location and read counts for all samples.
##' @author Jean Monlong
##' @export
quick.count <- function(files.df, bins.df, nb.cores=1,col.files=NULL,...){
  if(is.null(col.files)){
    if(nb.cores>1){
      bc.l = parallel::mclapply(files.df$sample, function(samp){
        bc.s = bin.bam(files.df$bam[files.df$sample==samp], bins.df, appendIndex.outfile=FALSE, ...)
        bc.s$bc$bc
      },mc.cores=nb.cores)
    } else {
      bc.l = lapply(files.df$sample, function(samp){
        bc.s = bin.bam(files.df$bam[files.df$sample==samp], bins.df, appendIndex.outfile=FALSE, ...)
        bc.s$bc$bc
      })
    }
    bc.df = createEmptyDF(c(sapply(c("chr","start","end"), function(cn)class(bins.df[,cn])), rep("integer",nrow(files.df))),nrow(bins.df))
  } else {
    if(nb.cores>1){
      bc.l = parallel::mclapply(files.df$sample, function(samp){
        bc.s = read.bedix(files.df[files.df$sample==samp, col.file], bins.df)
        bc.s[,4]
      },mc.cores=nb.cores)
    } else {
      bc.l = lapply(files.df$sample, function(samp){
        bc.s = read.bedix(files.df[files.df$sample==samp, col.file], bins.df)
        bc.s[,4]
      })
    }
    bc.df = createEmptyDF(c(sapply(c("chr","start","end"), function(cn)class(bins.df[,cn])), rep("numeric",nrow(files.df))),nrow(bins.df))
  }
  colnames(bc.df) = c("chr","start","end",as.character(files.df$sample))
  bc.df[,1:3] = bins.df[,c("chr","start","end")]
  for(ii in 1:length(bc.l)){
    bc.df[,3+ii] = bc.l[[ii]]
  }
  return(bc.df)
}
