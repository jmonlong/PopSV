##' The Z-score is computed by substracting the bin count by the average bin count
##' across the reference samples and dividing by their standard deviation. If
##' 'z.poisson' is TRUE, a score using Poisson distribution is also computed, using
##' the average bin count as an estimator of the lambda. Then the score with the lowest
##' absolute value is kept. This hybrid Z-score is to be used when some regions have low
##' coverage where it is more robust to use Poisson assumptions.
##' @title Single sample TMM normalization and test
##' @param test.sample the name of the sample to test.
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bc.gc.bg' are required and should be
##' present  after running 'initFileNames' function. Files should exist if
##' 'correct.GC' was run.
##' @param cont.sample the name of the sample used as control for the normalization.
##' @param bc.ref.f the path to the input file used for normalization (Optional).
##' @param norm.stats.f the name of the file with the statistic of the normalization run.
##' @param write.out.file should the result be written in files (from 'z' and 'fc' columns in 'files.df'). Default is TRUE.
##' @param compress.index should the output files be compressed and indexed. Default is TRUE.
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @param append should the results be appended to existing files. Default is FALSE.
##' @param col.file the column name in 'files.df' which inform which file contains the bin counts to use. Default is 'bc.gc.gz'.
##' @return a data.frame with columns :
##' \item{chr, start, end}{the location of the bin}
##' \item{bc}{the normalized bin count}
##' \item{z}{the Z-scores}
##' \item{fc}{the fold-change compared to the average bin count in the reference samples}
##' @author Jean Monlong
##' @export
tmm.test.sample <- function(test.sample, files.df, cont.sample, bc.ref.f=NULL, norm.stats.f, write.out.file=TRUE, compress.index=TRUE, z.poisson = FALSE, append=FALSE, col.file="bc.gc.gz") {
  chunk.size = 1000
  test.bc = utils::read.table(files.df[which(files.df$sample == test.sample),col.file], colClasses = c("character", "integer", "integer", "numeric"), header = TRUE)
  test.bins.order = paste(test.bc$chr, as.integer(test.bc$start), sep = "-")
  if(!is.null(bc.ref.f)){
    ref.headers = utils::read.table(bc.ref.f, nrows = 1, as.is = TRUE)
    colC = ifelse(ref.headers == cont.sample, "numeric", "NULL")
    colC[1:3] = c("character", "numeric", "numeric")
    cont.bc = utils::read.table(bc.ref.f, colClasses = colC, header = TRUE)
    colnames(cont.bc)[4] = "bc"
  } else {
    cont.bc = utils::read.table(files.df[which(files.df$sample == cont.sample), col.file], colClasses = c("character", "integer", "integer", "numeric"), header = TRUE)
  }
  id.cont = 1:nrow(cont.bc)
  names(id.cont) = paste(cont.bc$chr, as.integer(cont.bc$start), sep = "-")

  if (z.poisson) {
    z.comp.f <- function(x, mean.c, sd.c) {
      z.n = (x - mean.c)/sd.c
      z.p = stats::qnorm(stats::ppois(x, mean.c))
      n.ii = abs(z.n) < abs(z.p)
      z.p[which(n.ii)] = z.n[which(n.ii)]
      z.p
    }
  } else {
    z.comp.f <- function(x, mean.c, sd.c) {
      (x - mean.c)/sd.c
    }
  }

  ## prepare the output data.frame
  res.df = createEmptyDF(c("character", rep("integer", 2), rep("numeric", 3)),
    nrow(test.bc))
  colnames(res.df) = c("chr", "start", "end", "bc", "z", "fc")
  res.df$chr = test.bc$chr
  res.df$start = test.bc$start
  res.df$end = test.bc$end

  ## TMM normalize the test sample
  test.bc = cbind(test.bc, cont.bc$bc[id.cont[test.bins.order]])
  colnames(test.bc)[ncol(test.bc)] = cont.sample
  test.bc.normed = tmm.norm(test.bc, cont.sample, norm.stats.comp=FALSE)
  res.df$bc = test.bc.normed$bc.norm$bc
  rm(test.bc.normed)
  gc()

  ## Compute the Z-score and fold-change
  norm.stats.df = utils::read.table(norm.stats.f, as.is=TRUE, header=TRUE, sep='\t')
  id.norm.stats = 1:nrow(norm.stats.df)
  names(id.norm.stats) = paste(norm.stats.df[, 1], as.integer(norm.stats.df[, 2]), sep = "-")
  norm.stats.df = norm.stats.df[id.norm.stats[test.bins.order],]
  res.df$z = z.comp.f(res.df$bc, norm.stats.df$m, as.numeric(norm.stats.df$sd))
  res.df$fc = res.df$bc / norm.stats.df$m

  ## eventually write in files
  if(write.out.file){
    res.df = res.df[order(res.df$chr, res.df$start),]
    files.out = c(files.df[which(files.df$sample == test.sample), "z"],
                  files.df[which(files.df$sample == test.sample), "fc"],
                  files.df[which(files.df$sample == test.sample), "bc.gc.norm"])
    utils::write.table(res.df[,c("chr","start","end","z")], file = files.out[1], row.names = FALSE, quote = FALSE, sep = "\t", append=append, col.names=!append)
    utils::write.table(res.df[,c("chr","start","end","fc")], file = files.out[2], row.names = FALSE, quote = FALSE, sep = "\t", append=append, col.names=!append)
    utils::write.table(res.df[,c("chr","start","end","bc")], file = files.out[3], row.names = FALSE, quote = FALSE, sep = "\t", append=append, col.names=!append)
    if (compress.index) {
        comp.index.files(files.out)
    }
    res.df = files.out
  }

  return(res.df)
}
