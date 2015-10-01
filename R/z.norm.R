##' The Z-score is normalized to ensure centered distribution and similar variance across different GC classes. This is usually not necessary when the quality is comparable across the samples. However if some samples contain experiment-specific biases this approach can reduce false calls. 
##' @title Z-score normalization
##' @param z.df a data.frame with the Z-scores. 
##' @param gc.df a data.frame with the GC content.
##' @param class.size the number of bins in each GC classes. Default is 5000.
##' @return a vector with the normalized Z-scores.
##' @author Jean Monlong
z.norm <- function(z.df,gc.df, class.size=5000){
  ## Slightly faster merge
  mergeGC <- function(z.df, bins.df){
    gcc = bins.df$GCcontent
    names(gcc) = with(bins.df, paste(chr,as.integer(start),as.integer(end)))
    z.df$GCcontent = with(z.df, gcc[paste(chr,as.integer(start),as.integer(end))])
    z.df
  }

  z.df = mergeGC(z.df, gc.df)
  z.df$gc.class = cut(z.df$GCcontent, include.lowest = TRUE, breaks=unique(quantile(z.df$GCcontent,probs=seq(0,1,class.size/nrow(z.df)))))
  msd.l = tapply(z.df$z, z.df$gc.class, function(z)c(median(z, na.rm=TRUE), mad(z, na.rm=TRUE)))
  m.v = unlist(lapply(msd.l, "[", 1))
  sd.v = unlist(lapply(msd.l, "[", 2))
  (z.df$z-m.v[as.character(z.df$gc.class)])/sd.v[as.character(z.df$gc.class)]
}
