##' Split the bins into big and small chunks. A big chunk represents all the bins used to
##' normalize bins of this chunk. Small chunk represents the bins to analyze by one job
##' on the cluster. While the size of the small chunks is not important and can be adjusted
##' to fit the cluster, the big chunk size will impact on the efficiency of the normalization
##' (the bigger the better).
##' @title Split the bins in chunks for parallel normalization
##' @return an updated data.frame with new columns 'sm.chunk', 'bg.chunk' and 'bin' with the small chunk
##' ID, big chunk ID and bin definition.
##' @author Jean Monlong
##' @param bins.df a data.frame with the bins definition (one row per bin). E.g. created by
##' 'fragment.genome.hg19'.
##' @param bg.chunk.size the number of bins in a big chunk.
##' @param sm.chunk.size the number of bins in a small chunk.
##' @param large.chr.chunks should the big chunks be made of a few large genomic sub-regions ? Default is false. Normalization is faster (but a bit less efficient) than when using random bins. Recommended when dealing with a large number of bins.
##' @export
chunk.bin <- function(bins.df, bg.chunk.size = 1e+05, sm.chunk.size = 1000, large.chr.chunks = FALSE) {
  bins.df = bins.df[order(as.character(bins.df$chr), bins.df$start),]
  if (large.chr.chunks) {
    nb.supchunks = 50
    bins.df$bg.chunk = rep(sample(rep(1:ceiling(nrow(bins.df)/bg.chunk.size), nb.supchunks)), each = ceiling(bg.chunk.size/nb.supchunks))[1:nrow(bins.df)]
  } else {
    bins.df$bg.chunk = sample(rep(1:ceiling(nrow(bins.df)/bg.chunk.size), each = bg.chunk.size)[1:nrow(bins.df)])
  }
  bg.chunk = chr = NULL  ## Uglily appease R checks
  bins.df = dplyr::mutate(dplyr::group_by(bins.df, bg.chunk),
                          sm.chunk = paste(bg.chunk,
                                           sample(rep(1:ceiling(length(chr)/sm.chunk.size),
                                                      each = sm.chunk.size)[1:length(chr)]),
                                           sep = "-"))
  bins.df$bin = paste(bins.df$chr, bins.df$start, sep = "-")
  bins.df = dplyr::ungroup(bins.df)
  bins.df$start = as.integer(bins.df$start)
  bins.df$end = as.integer(bins.df$end)
  return(bins.df)
}
