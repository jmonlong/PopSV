##' Bin counts are normalized one bin at a time, using a subset of the bins that look
##' similar across the reference samples.
##'
##' @title Targeted-normalization of bin counts
##' @param bc a matrix or data.frame with the bin counts (bin x sample).
##' @param cont.sample the sample to use as baseline for the pairwise normalization.
##' All the samples will be normalized to it.
##' @param nb.support.bins the number of bins to use for the normalization.
##' @param bins a vector the names of the bins to normalize. If NULL (default), all
##' bins are normalized.
##' @param save.support.bins if TRUE (default) the bins used for the normalization are
##' saved in the output object 'norm.stats'.
##' @param bootstrap should the supporting bins defined using a bootstrap approach. Default is FALSE.
##' @return a list with
##' \item{norm.stats}{a data.frame witht some metrics about the normalization of each
##' bin (row) : correlation with worst supporting bin; coverage average and standard
##' deviation; number of outlier reference samples; supporting bins.}
##' \item{bc.norm}{a matrix with the normalized bin counts (bin x sample).}
##' \item{nb.support.bins, cont.sample, z.poisson}{a backup of the input parameters.}
##' \item{cont.sample}{the name of the sample used to normalize all samples.}
##' @author Jean Monlong
##' @export
tn.norm <- function(bc, cont.sample, nb.support.bins = 1000, bins = NULL, save.support.bins = TRUE, bootstrap = FALSE) {
  
  all.samples = setdiff(colnames(bc), c("chr", "start", "end"))
  rownames(bc) = paste(bc$chr, as.integer(bc$start), sep = "-")
  if (is.null(bins)) 
    bins = rownames(bc)
  
  if (save.support.bins) {
    norm.stats = createEmptyDF(c("character", rep("integer", 2), rep("numeric", 
      4), rep("character", nb.support.bins)), length(bins))
    colnames(norm.stats) = c("chr", "start", "end", "m", "sd", "nb.remove", "d.max", 
              paste("b", 1:nb.support.bins, sep = ""))
  } else {
    norm.stats = createEmptyDF(c("character", rep("integer", 2), rep("numeric", 
      4)), length(bins))
    colnames(norm.stats) = c("chr", "start", "end", "m", "sd", "nb.remove", "d.max")
  }
  bc.norm = createEmptyDF(c("character", rep("integer", 2), rep("numeric", length(all.samples))), 
    length(bins))
  colnames(bc.norm) = c("chr", "start", "end", all.samples)
  norm.stats$chr = bc.norm$chr = bc[bins, "chr"]
  norm.stats$start = bc.norm$start = bc[bins, "start"]
  norm.stats$end = bc.norm$end = bc[bins, "end"]
  
  bc = t(as.matrix(bc[, all.samples]))
  for (bin.ii in 1:length(bins)) {
    
    bin = bins[bin.ii]
    bc.i = bc[, bin]
    if (any(!is.na(bc.i) & bc.i != 0)) {
      if (sum(as.numeric(bc.i) > 0) < 5) {
        d.o.i = sample(1:ncol(bc), nb.support.bins)
        d.max = -1
      } else {
        if(bootstrap){
          d.i = sapply(1:5,function(ii){
            bs.samps = sample(1:length(bc.i), length(bc.i)*.5)
            1 - as.numeric(suppressWarnings(cor(as.numeric(bc.i[bs.samps]), bc[bs.samps,], use = "pairwise.complete.obs")))
          })
          d.i = apply(d.i, 1, function(ee) suppressWarnings(max(ee, na.rm=TRUE)))
          if(any(is.infinite(d.i))) {
            d.i[which(is.infinite(d.i))] = NA
          }
        } else {
          d.i = 1 - as.numeric(suppressWarnings(cor(as.numeric(bc.i), bc, use = "pairwise.complete.obs")))
        }
        d.o.i = order(d.i)[1:nb.support.bins]
        d.max = as.numeric(d.i[d.o.i[nb.support.bins]])
      }
      if (!is.infinite(d.max) & !is.na(d.max)) {
        bc.g = t(bc[, d.o.i])
        bin.for.norm = colnames(bc)[d.o.i]
        norm.coeff = norm.tm.opt(bc.g, ref.col = bc.g[, cont.sample])
        bc.t = bc.i * norm.coeff
        msd = mean.sd.outlierR(bc.t, 1e-06)
        if (any(!is.na(bc.t))) {
          norm.stats[bin.ii, 4:7] = round(c(msd$m, msd$sd, msd$nb.remove, d.max), 3)
          if (save.support.bins) 
            norm.stats[bin.ii, 8:ncol(norm.stats)] = bin.for.norm
          bc.norm[bin.ii, -(1:3)] = round(bc.t, 2)
        }
      }
    }
    
  }
  
  return(list(norm.stats = norm.stats, bc.norm = bc.norm, nb.support.bins = nb.support.bins, 
              cont.sample = cont.sample))
} 
