##' Fragment hg19 genome into consecutive bins of equal size. This function
##' require package \code{BSgenome.Hsapiens.UCSC.hg19} to be installed.
##' @title Fragment hg19 genome
##' @param bin.size the size of the bins
##' @return a data.frame with columns 'chr', 'start' and 'end'.
##' @author Jean Monlong
##' @export
fragment.genome.hp19 <- function(bin.size = 1000) {
    if (!require(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)) {
        stop("Please install BSgenome first by running:\n> source(\"http://bioconductor.org/biocLite.R\")\n> biocLite(\"BSgenome.Hsapiens.UCSC.hg19\")\n")
    }
    seql.1.22 = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[paste0("chr", 
        1:22)]
    fragment.chr <- function(chr.i) {
        starts = as.integer(seq(1, seql.1.22[chr.i], bin.size))
        ends = as.integer(starts + bin.size - 1)
        ends[length(ends)] = as.integer(seql.1.22[chr.i])
        data.frame(chr = chr.i, start = starts, end = ends)
    }
    plyr::ldply(lapply(1:22, fragment.chr), identity)
} 
