##' Quantile normalization of the bin counts. After such normalization, each sample has exactly the same global bin count distribution.
##' @title Quantile normalization
##' @param bc.df a data.frame with 'chr', 'start', 'end' columns and then one column per sample with its bin counts.
##' @param nb.cores the number of cores to use. Default is 1.
##' @param norm.stats.comp Should some statistics on the normalized bin count be computed (mean, sd, outliers). Default is TRUE.
##' @return a list with
##' \item{norm.stats}{a data.frame witht some metrics about the normalization of each
##' bin (row) : coverage average and standard deviation; number of outlier samples}
##' \item{bc.norm}{a data.frame, similar to the input 'bc.df', with the normalized bin counts.}
##' @author Jean Monlong
##' @export
quant.norm <- function(bc.df, nb.cores = 1, norm.stats.comp = TRUE) {
  all.samples = setdiff(colnames(bc.df), c("chr", "start", "end"))
  rownames(bc.df) = bins = paste(bc.df$chr, as.integer(bc.df$start), as.integer(bc.df$end), sep = "-")

  bc.norm = createEmptyDF(c("character", rep("integer", 2), rep("numeric", length(all.samples))), length(bins))
  colnames(bc.norm) = c("chr", "start", "end", all.samples)
  bc.norm$chr = bc.df$chr
  bc.norm$start = bc.df$start
  bc.norm$end = bc.df$end

  bc.mat = as.matrix(bc.df[, all.samples])
  bc.cov = as.numeric(parallel::mclapply(1:ncol(bc.mat), function(cc) stats::median(bc.mat[, cc], na.rm = TRUE), mc.cores = nb.cores))
  bc.mat = (bc.mat * stats::median(bc.cov)) %*% diag(1/bc.cov)
  bc.sorted = matrix(as.numeric(unlist(parallel::mclapply(sample(1:ncol(bc.mat),5), function(cc) sort(bc.mat[, cc]), mc.cores = nb.cores))), nrow(bc.mat))
  med = round(as.numeric(parallel::mclapply(1:nrow(bc.sorted), function(rr) stats::median(bc.sorted[rr, ], na.rm = TRUE), mc.cores = nb.cores)), 2)
  bc.norm[, -(1:3)] = matrix(as.numeric(unlist(parallel::mclapply(1:ncol(bc.mat),
           function(cc) med[rank(bc.mat[, cc])], mc.cores = nb.cores))), nrow(bc.mat))

  if (norm.stats.comp) {
    norm.stats = createEmptyDF(c("character", rep("integer", 2), rep("numeric", 3)), length(bins))
    colnames(norm.stats) = c("chr", "start", "end", "m", "sd", "nb.remove")
    norm.stats$chr = bc.df$chr
    norm.stats$start = bc.df$start
    norm.stats$end = bc.df$end
    norm.stats[, 4:6] = matrix(as.numeric(unlist(parallel::mclapply(1:nrow(bc.norm),
                function(rr) {
                  msd = mean.sd.outlierR(as.numeric(bc.norm[rr, all.samples]))
                  return(c(msd$m, msd$sd, msd$nb.remove))
                }, mc.cores = nb.cores))), nrow(bc.norm))
  } else {
    norm.stats = NULL
  }

  return(list(norm.stats = norm.stats, bc.norm = bc.norm))
}
