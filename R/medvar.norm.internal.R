##' The median and deviation to the median is normalized across samples. 
##' @title Median-variance normalization of bin counts
##' @param bc a matrix or data.frame with the bin counts (bin x sample).
##' @return a matrix with the normalized bin counts (bin x sample).
##' @author Jean Monlong
##' @keywords internal
medvar.norm.internal <- function(bc) {
  med = apply(bc, 2, median, na.rm = TRUE)
  if(any(med==0)) med[which(med==0)] = 1
  med.c = mean(med)
  bc = t(t(bc) * med.c/med)
  bc = bc - med.c
  md = apply(bc, 2, function(x) median(abs(x), na.rm = TRUE))
  md.c = median(abs(bc), na.rm = TRUE)
  bc = t(t(bc) * md.c/md)
  return(bc + med.c)
} 
