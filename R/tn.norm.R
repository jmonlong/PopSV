##' Bin counts are normalized one bin at a time, using a subset of the bins that look
##' similar across the reference samples.
##'
##' A specific set of bins, defined by 'bins=', can de normalized using all the bins in 'bc'. The bin names in 'bins=' should be the chr and start position as in '1-501' (for chr 1, start 501).
##'
##' The default approach ('norm="1pass"') looks for supporting across all samples. A more robust approach, finds supporting bins after trimming a few outlier samples (potential CNV). 'norm="trim"' normalization is better and recommended for small bins.
##' @title Targeted-normalization of bin counts
##' @param bc a matrix or data.frame with the bin counts (bin x sample).
##' @param cont.sample the sample to use as baseline for the pairwise normalization.
##' All the samples will be normalized to it.
##' @param nb.support.bins the number of bins to use for the normalization.
##' @param bins a vector the names of the bins to normalize. If NULL (default), all
##' bins are normalized.
##' @param save.support.bins if TRUE (default) the bins used for the normalization are
##' saved in the output object 'norm.stats'.
##' @param norm the type of normalization. '1pass' (default) means one pass of normalization. Other options is 'trim'.
##' @param force.diff.chr should the supporting bins be forced to be in a different chromosome. Default is TRUE.
##' @return a list with
##' \item{norm.stats}{a data.frame witht some metrics about the normalization of each
##' bin (row) : correlation with worst supporting bin ('d.max'); coverage average ('m') and standard deviation ('sd'); number of outlier reference samples ('nb.remove'); supporting bins.}
##' \item{bc.norm}{a matrix with the normalized bin counts (bin x sample).}
##' \item{nb.support.bins, cont.sample}{a backup of the input parameters.}
##' @author Jean Monlong
##' @export
tn.norm <- function(bc, cont.sample, nb.support.bins = 1000, bins = NULL, save.support.bins = TRUE, norm = c("1pass", "trim"), force.diff.chr=TRUE) {

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
  bc.norm = createEmptyDF(c("character", rep("integer", 2), rep("numeric", length(all.samples))), length(bins))
  colnames(bc.norm) = c("chr", "start", "end", all.samples)
  norm.stats$chr = bc.norm$chr = bc[bins, "chr"]
  norm.stats$start = bc.norm$start = bc[bins, "start"]
  norm.stats$end = bc.norm$end = bc[bins, "end"]

  chrs = bc$chr
  denorm.factor = runif(length(all.samples), 1,1.5)
  denorm.factor[which(all.samples==cont.sample)] = 1
  bc = t(as.matrix(bc[, all.samples]))
  bc = denorm.factor * bc

  trimmed.index <- function(x,nb.trim=5){
    xs = sort(x)
    which(x>xs[nb.trim] & x<xs[length(x)-nb.trim+1])
  }

  for (bin.ii in 1:length(bins)) {
    bin = bins[bin.ii]
    bc.i = bc[, bin]
    if (any(!is.na(bc.i) & bc.i != 0)) {
      if (sum(as.numeric(bc.i) > 0) < 5) {
        d.o.i = sample(1:ncol(bc), nb.support.bins)
        d.max = -1
      } else {
        if(norm[1]=="trim"){
          trim.i = trimmed.index(bc.i/denorm.factor, 3)
          d.i = 1 - as.numeric(suppressWarnings(cor(as.numeric(bc.i[trim.i]), bc[trim.i,], method="spearman", use = "pairwise.complete.obs")))
        } else {
          d.i = 1 - as.numeric(suppressWarnings(cor(as.numeric(bc.i), bc, use = "pairwise.complete.obs")))
        }
        if(force.diff.chr){
          chr.i = bc.norm$chr[bin.ii]
          if(any(chr.i==chrs)){
            d.i[which(chr.i==chrs)] = NA
          }
        }
        d.o.i = order(d.i)[1:nb.support.bins]
        d.max = as.numeric(d.i[d.o.i[nb.support.bins]])
      }
      if (!is.infinite(d.max) & !is.na(d.max)) {
        bc.g = t(bc[, d.o.i])
        bin.for.norm = colnames(bc)[d.o.i]
        norm.coeff = norm.tm.opt(bc.g, ref.col = bc.g[, cont.sample])
        bc.t = bc.i * norm.coeff
        msd = mean.sd.outlierR(bc.t)
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
