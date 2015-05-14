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
    bc.l = parallel::mclapply(files.df$sample, function(samp) {
      bc.s = bin.bam(files.df$bam[files.df$sample == samp], bins.df, appendIndex.outfile = FALSE, 
        ...)
      bc.s$bc$bc
    }, mc.cores = nb.cores)
  } else {
    bc.l = parallel::mclapply(files.df$sample, function(samp) {
      bc.s = read.bedix(files.df[files.df$sample == samp, col.files], bins.df,exact.match=TRUE)
      bc.s[, 4]
    }, mc.cores = nb.cores)
  }
  bc.df = matrix(unlist(bc.l), ncol=length(bc.l))
  colnames(bc.df) =  as.character(files.df$sample)
  return(cbind( bins.df[, c("chr", "start", "end")], bc.df))
} 
