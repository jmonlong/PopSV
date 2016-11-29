##' Find bins that are not covered and remove them from the input 'bins.df'. This function is useful for targeted sequencing when the binning was unaware of the targeted regions. In practice, the bin counts across a few samples are used to assess if a bin is covered or not.
##'
##' If 'bc.med.min' is NULL, the function will quickly model the coverage due to off-target reads and set the minimum threshold as (median+3sd) of the off-target coverage. The off-tagert coverage is approximated by the maximum non-null coverage in bins with median coverage of 0. If a bin has 0 read in more than half of the samples, it's likely a non-covered bins. However sometimes a few samples have some reads there. The distribution of these noise-coverage is used to fix the minimum threshold in median coverage.
##' @title Filter non-covered bins
##' @param bins.df the data.frame with the bins information that will be subsetted.
##' @param files.df the data.frame with the file information, e.g. location of the bin count files.
##' @param nb.samples the number of samples to use. Default is 10.
##' @param bc.med.min If non-NULL, the minimum median coverage to use for the filtering. If NULL (default), this value will be estimated.
##' @param plot Display a graph ? Default is TRUE.
##' @return a subset of the input 'bins.df' data.frame
##' @author Jean Monlong
##' @export
filter.noncovered.bins <- function(bins.df, files.df, nb.samples=10, bc.med.min=NULL, plot=TRUE){
  bc.l = lapply(sample.int(nrow(files.df),nb.samples), function(ii){
    bc = read.table(files.df$bc.gz[ii], as.is=TRUE, header=TRUE)
    bc$sample = files.df$sample[ii]
    bc
  })

  bc.mat = matrix(unlist(lapply(bc.l, function(df)df$bc)), nrow(bc.l[[1]]))
  colnames(bc.mat) = unlist(lapply(bc.l, function(df)df$sample[1]))
  
  bc.s = data.frame(bc.med = apply(bc.mat, 1, median), bc.max = apply(bc.mat, 1, max))
  
  if(is.null(bc.med.min)){
    bc.med.min = median(bc.s$bc.max[which(bc.s$bc.med==0 & bc.s$bc.max!=0)])
    bc.med.min.sd = mad(bc.s$bc.max[which(bc.s$bc.med==0 & bc.s$bc.max!=0)])
    bc.med.min = bc.med.min + 3*bc.med.min.sd
  }
  
  if(plot){
    winsor <- function(x,u=100){
      if(any(x>u)) x[x>u] = u
      x
    }
    ## Uglily appeases R checks
    bc.max.w = NULL
    bc.s$bc.max.w = winsor(bc.s$bc.max,2*bc.med.min)
    print(ggplot2::ggplot(bc.s[which(bc.s$bc.med==0),], ggplot2::aes(x=bc.max.w)) + ggplot2::geom_histogram(binwidth=1) + ggplot2::theme_bw() + ggplot2::xlab("maximum bin count") + ggplot2::ggtitle("Bins with 0 median bin count") + ggplot2::ylab("bin") + ggplot2::geom_vline(xintercept=bc.med.min, linetype=2))
  }
  
  bins = paste0(bc.l[[1]]$chr, "-", bc.l[[1]]$start)
  ## Check that the order is the same in the files
  check.order = lapply(bc.l, function(df) all(bins==paste0(df$chr, "-", df$start)))
  if(!all(unlist(check.order))){
    stop("Different bin order in the bin count files. This won't work...")
  }
  
  covered.bins = bins[which(bc.s$bc.med >= bc.med.min)]
  covered.bins = which(paste0(bins.df$chr, "-", bins.df$start) %in% covered.bins)
  bins.df[covered.bins,]
}
