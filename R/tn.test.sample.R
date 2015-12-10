##' Bin counts from one sample are normalized following instructions from
##' a previous targeted-normalization run.
##'
##' The Z-score is computed by substracting the bin count by the average bin count
##' across the reference samples and dividing by their standard deviation. If
##' 'z.poisson' is TRUE, a score using Poisson distribution is also computed, using
##' the average bin count as an estimator of the lambda. Then the score with the lowest
##' absolute value is kept. This hybrid Z-score is to be used when some regions have low
##' coverage where it is more robust to use Poisson assumptions.
##' @title Single sample targeted normalization and test
##' @param test.sample the name of the sample to test.
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bc.gc.bg' are required and should be
##' present  after running 'initFileNames' function. Files should exist if
##' 'correct.GC' was run.
##' @param cont.sample the name of the sample used as control for the normalization.
##' @param bc.ref.f the path to the input file used for targeted normalization ('tn.norm').
##' @param norm.stats.f the name of the file with the statistic of the targeted normalization run.
##' @param write.out.file should the result be written in files (from 'z' and 'fc' columns in 'files.df'). Default is TRUE.
##' @param compress.index should the output files be compressed and indexed. Default is TRUE.
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @param aberrant.cases if TRUE (default) a more robust (but sligthly longer) normalization
##' is performed on cases to deal with potential large chromosomal aberrations. In practice,
##' it is recommended for cancer but can be turned off if less than ~20% of the genome is expected
##' to be affected.
##' @param append should the results be appended to existing files. Default is FALSE.
##' @param col.file the column name in 'files.df' which inform which file contains the bin counts to use. Default is 'bc.gc.gz'.
##' @return a data.frame with columns :
##' \item{chr, start, end}{the location of the bin}
##' \item{bc}{the normalized bin count}
##' \item{z}{the Z-scores}
##' \item{fc}{the fold-change compared to the average bin count in the reference samples}
##' @author Jean Monlong
##' @export
tn.test.sample <- function(test.sample, files.df, cont.sample, bc.ref.f=NULL, norm.stats.f, write.out.file=TRUE, compress.index=TRUE, z.poisson = FALSE, aberrant.cases = FALSE, append=FALSE, col.file="bc.gc.gz") {
  chunk.size = 1000
  test.bc = read.table(files.df[which(files.df$sample == test.sample),col.file], colClasses = c("character", "integer", "integer", "numeric"), header = TRUE)
  id.test = 1:nrow(test.bc)
  names(id.test) = paste(test.bc$chr, as.integer(test.bc$start), sep = "-")
  if(!is.null(bc.ref.f)){
    ref.headers = read.table(bc.ref.f, nrows = 1, as.is = TRUE)
    colC = ifelse(ref.headers == cont.sample, "numeric", "NULL")
    colC[1:3] = c("character", "integer", "integer")
    cont.bc = read.table(bc.ref.f, colClasses = colC, header = TRUE)
    colnames(cont.bc)[4] = "bc"
  } else {
    cont.bc = read.table(files.df[which(files.df$sample == cont.sample), col.file], colClasses = c("character", "integer", "integer", "numeric"), header = TRUE)
  }
  id.cont = 1:nrow(cont.bc)
  names(id.cont) = paste(cont.bc$chr, as.integer(cont.bc$start), sep = "-")

  if (z.poisson) {
    z.comp.f <- function(x, mean.c, sd.c) {
      z.n = (x - mean.c)/sd.c
      z.p = qnorm(ppois(x, mean.c))
      n.ii = abs(z.n) < abs(z.p)
      z.p[which(n.ii)] = z.n[which(n.ii)]
      z.p
    }
  } else {
    z.comp.f <- function(x, mean.c, sd.c) {
      (x - mean.c)/sd.c
    }
  }

  res.df = createEmptyDF(c("character", rep("integer", 2), rep("numeric", 3)),
    nrow(test.bc))
  colnames(res.df) = c("chr", "start", "end", "bc", "z", "fc")
  res.df$chr = test.bc$chr
  res.df$start = test.bc$start
  res.df$end = test.bc$end

  if (aberrant.cases) {
    norm.bin <- function(bc.bin.arg, test.bc.arg, cont.bc.arg, bc.mean.arg=NULL, chrs.arg=NULL) {
      norm.coeff = norm.tm.opt(as.matrix(test.bc.arg),
        ref.col = cont.bc.arg, bc.mean.norm = bc.mean.arg,
        chrs = chrs.arg)
      bc.bin.arg * norm.coeff
    }
  } else {
    norm.bin <- function(bc.bin.arg, test.bc.arg, cont.bc.arg, bc.mean.arg=NULL, chrs.arg=NULL) {
      norm.coeff = norm.tm.opt(as.matrix(test.bc.arg),
        ref.col = cont.bc.arg)
      bc.bin.arg * norm.coeff
    }
  }

  con = file(norm.stats.f, "r")
  headers = unlist(strsplit(readLines(con, n = 1), "\t"))
  while (length((lines = readLines(con, n = chunk.size))) > 0) {
    norm.chunk = matrix(unlist(strsplit(lines, "\t")), ncol = length(headers),
      byrow = TRUE)
    colnames(norm.chunk) = headers
    norm.chunk = norm.chunk[which(norm.chunk[, 4] != "NA"), ]

    test.bc.chunk = matrix(test.bc$bc[id.test[unlist(norm.chunk[,-(1:7)])]],nrow(norm.chunk))
    cont.bc.chunk = matrix(cont.bc$bc[id.cont[unlist(norm.chunk[,-(1:7)])]],nrow(norm.chunk))

    chrs.chunk=NULL
    if(aberrant.cases){
      chrs.chunk = matrix(test.bc$chr[id.test[unlist(norm.chunk[,-(1:7)])]],nrow(norm.chunk))
    }

    id.bins = id.test[paste(norm.chunk[, 1], as.integer(norm.chunk[, 2]), sep = "-")]
    bc.mean.chunk = as.numeric(norm.chunk[,4])
    res.df$bc[id.bins] = sapply(1:nrow(norm.chunk), function(ii){
      norm.bin(test.bc$bc[id.bins[ii]], test.bc.chunk[ii,], cont.bc.chunk[ii,], bc.mean=bc.mean.chunk[ii], chrs=chrs.chunk[ii,])
    })
    res.df$z[id.bins] = z.comp.f(res.df$bc[id.bins], bc.mean.chunk, as.numeric(norm.chunk[,5]))
    res.df$fc[id.bins] = res.df$bc[id.bins] / bc.mean.chunk
  }
  close(con)

  if(write.out.file){
    files.out = c(files.df[which(files.df$sample == test.sample), "z"], files.df[which(files.df$sample == test.sample), "fc"], files.df[which(files.df$sample == test.sample), "bc.gc.norm"])
    write.table(res.df[,c("chr","start","end","z")], file = files.out[1], row.names = FALSE, quote = FALSE, sep = "\t", append=append, col.names=!append)
    write.table(res.df[,c("chr","start","end","fc")], file = files.out[2], row.names = FALSE, quote = FALSE, sep = "\t", append=append, col.names=!append)
    write.table(res.df[,c("chr","start","end","bc")], file = files.out[3], row.names = FALSE, quote = FALSE, sep = "\t", append=append, col.names=!append)
    if (compress.index) {
      comp.index.files(files.out)
    }
    res.df = files.out
  }

  return(res.df)
}
