##' Finds the local maximum.
##' @title Local maximum
##' @param x a vector with numeric values
##' @param min.max.prop the minimum height of the maximum to be considered. Default is 0.1, i.e. all maximum at least 10% as height as the highest maximum.
##' @return a list with the local maximum position and height.
##' @author Jean Monlong
##' @keywords internal
localMax <- function(x, min.max.prop = 0.1) {
  d = density(x, na.rm = TRUE)
  im = 1 + which(diff(sign(diff(d$y))) == -2)
  my = max(d$y)
  max.id = im[which(d$y[im] >= min.max.prop * my)]
  max.id.o = max.id[order(d$y[max.id], decreasing = TRUE)]
  return(list(lM = d$x[max.id], h = d$y[max.id]/my))
}
