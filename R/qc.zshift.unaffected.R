##' Checks that the Z-scores of unaffected samples are centered in 0. When controls are not perfectly matched, the Z-scores of cases can be slightly shifted and produce false calls. This function computes a measure of such bias. 
##' @title QC: Z-score bias in unaffected samples
##' @param files.df a data.frame with the path to relevant files for each samples.
##' @param bins.df a data.frame with information on which bins to analyze.
##' @param nb.cores the number of cores to use. Default is 1.
##' @return a data.frame with a shift measure for each bin.
##' @author Jean Monlong
##' @export
qc.zshift.unaffected <- function(files.df, bins.df, nb.cores=1){
  z.df = quick.count(files.df, bins.df, nb.cores=nb.cores, col.files="z.gz")

  zshift <- function(z){
    lm = localMax(z)
    min(abs(lm$lM))
  }

  zs.df = z.df[,c("chr","start","end")]
  zs.df$zs = apply(z.df[,as.character(files.df$sample)], 1, zshift)

  return(zs.df)
}
