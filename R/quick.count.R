##' This function retrieves bin count for all samples in the input file information
##' data.frame. It was designed to count reads in a small amount of bins (a subset
##' of what will be used for the full analysis). This counts can be used later to
##' compare a potentially large number samples and guide the selection of samples
##' to consider in the analysis.
##'
##' If `col.files="bc.gc.gz"`, for example, the bin counts will be merged from the pre-computed files (the ones created by functions 'bin.bam' and 'correct.GC').
##' @title Counts reads across samples in a small number of bins
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bam' are required.
##' @param bins.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param nb.cores number of cores to use. If higher than 1, \code{parallel}
##' package is used to parallelized the counting.
##' @param col.files if non-NULL, the name of the column in 'files.df' with the filenames where the information should be retrieved from.
##' @param nb.rand.bins if non-NULL, the number of bins to randomly choose. Default is NULL.
##' @param ... extra paramaters for 'bin.bam' function
##' @return a data.frame with the bin location and read counts for all samples.
##' @author Jean Monlong
##' @export
quick.count <- function(files.df, bins.df, nb.cores = 1, col.files = NULL, nb.rand.bins=NULL, ...) {
  if(!is.null(nb.rand.bins)){
    bins.df = bins.df[sample.int(nrow(bins.df), min(nrow(bins.df), nb.rand.bins)),]
  }

  if (is.null(col.files)) {
    bins.df = bin.bam(files.df$bam[1], bins.df, appendIndex.outfile = FALSE, ...)$bc
    bc.l = parallel::mclapply(files.df$sample[-1], function(samp) {
      bc.s = bin.bam(files.df$bam[files.df$sample == samp], bins.df, appendIndex.outfile = FALSE, 
        ...)
      bc.s$bc$bc
    }, mc.cores = nb.cores)
    bc.l = c(list(bins.df$bc), bc.l)
  } else {
    bins.df = read.bedix(files.df[1, col.files], bins.df)
    bc.l = parallel::mclapply(files.df$sample[-1], function(samp) {
      bc.s = read.bedix(files.df[files.df$sample == samp, col.files], bins.df)
      bc.s[, 4]
    }, mc.cores = nb.cores)
    bc.l = c(list(bins.df[,4]), bc.l)
  }
  if(any(unlist(lapply(bc.l, length))!=nrow(bins.df))){
    stop("Inconsistent merging.", bc.l[[which(unlist(lapply(bc.l, length))!=nrow(bins.df))[1]]][1])
  }
  bc.df = matrix(unlist(bc.l), ncol=length(bc.l))
  colnames(bc.df) =  as.character(files.df$sample)
  return(cbind( bins.df[, c("chr", "start", "end")], bc.df))
} 
