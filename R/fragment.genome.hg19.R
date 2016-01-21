##' Fragment hg19 genome into consecutive bins of equal size. This function
##' require package \code{BSgenome.Hsapiens.UCSC.hg19} to be installed.
##' @title Fragment hg19 genome
##' @param bin.size the size of the bins
##' @param slid.window the size of the sliding window. Default is the size of the bin, i.e. not overlapping windows.
##' @param chr.prefix should chromosome name be in the form "chr1" instead of "1". Default is FALSE.
##' @param XY.chr should chromosome X and Y be included. Default is FALSE. If TRUE, male and female should be analyzed separately, at least for chromosomes X and Y.
##' @return a data.frame with columns 'chr', 'start' and 'end'.
##' @author Jean Monlong
##' @import GenomicRanges
##' @import GenomeInfoDb
##' @export
fragment.genome.hg19 <- function(bin.size = 1000, slid.window=bin.size, chr.prefix=FALSE, XY.chr=FALSE) {
  if (!require(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)) {
    stop("Please install BSgenome first by running:\n> source(\"http://bioconductor.org/biocLite.R\")\n> biocLite(\"BSgenome.Hsapiens.UCSC.hg19\")\n")
  }

  if(XY.chr){
    chrs = c(1:22, "X", "Y")
    message("Chromosome X/Y included : male and female should be analyzed separately.")
  } else {
    chrs = 1:22
  }

  chr.chrs = paste0("chr",chrs)
  if(chr.prefix){
    chrs = chr.chrs
  }

  seql.chrs = seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[chr.chrs]

  fragment.chr <- function(chr.i) {
    starts = as.integer(seq(1, seql.chrs[chr.i], slid.window))
    ends = as.integer(starts + bin.size - 1)
    if(any(ends>as.integer(seql.chrs[chr.i]))){
      ends[which(ends>as.integer(seql.chrs[chr.i]))] = as.integer(seql.chrs[chr.i])
    }
    data.frame(chr = chrs[chr.i], start = starts, end = ends)
  }

  do.call(rbind, lapply(1:length(chrs), fragment.chr))
}
