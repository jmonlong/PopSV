##' Chromosomes with clear aneuploidy are detected using a simple threshold-based approach. It is used to flag chromosome with complete aneuploidy that could handicap later analysis. It might not detect partial chromosomal aberration or aneuploidy in samples with noisy read coverage.
##' @title Flag chromosomal aneuploidy
##' @param samp the name of the sample to analyze.
##' @param files.df a data.frame with information about samples such as the path to the count files.
##' @param col.file the name of the column in 'files.df' with the path information to use.
##' @param nb.bins the number of bins to use when subsampling each chromosome. Default is 1000.
##' @param prop.aneu the proportion of normal bins under which a chromosome is considered aneuploid.
##' @return a vector with the names of the flagged(aneuploid) chromosomes.
##' @author Jean Monlong
##' @export
aneuploidy.flag <- function(samp, files.df, col.file="bc.bg", nb.bins=1e3, prop.aneu=.1){
  
  localMax <- function(x,min.max.prop=.1,max=FALSE){
    d = density(x,na.rm=TRUE)
    im = 1+which(diff(sign(diff(d$y)))==-2)
    my = max(d$y)
    max.id = im[which(d$y[im] >= min.max.prop * my)]
    max.id.o = max.id[order(d$y[max.id],decreasing=TRUE)]
    return(list(lM=d$x[max.id], h=d$y[max.id]/my))
  }

  df = dplyr::do(dplyr::group_by(subset(df, bc>0), chr), {.[sample.int(nb.bins),]})
  lm.o = localMax(df$bc)
  lm.o = lm.o$lM[which.max(lm.o$h)]
  core.chrs = df$chr[order(abs(lm.o-df$bc))[1:(nrow(df)*.3)]]
  cchrs.t = table(core.chrs)
  aneu.chrs = names(cchrs.t)[which(cchrs.t > prop.aneu*nb.bins*.3)]

  return(aneu.chrs)  
}
