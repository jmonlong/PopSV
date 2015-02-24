##' Fragment hg19 genome into consecutive bins of equal size. This function
##' require package \code{BSgenome.Hsapiens.UCSC.hg19} to be installed.
##' @title Fragment hg19 genome
##' @param bin.size the size of the bins
##' @return a data.frame with columns 'chr', 'start' and 'end'.
##' @author Jean Monlong
##' @export
fragment.genome.hp19 <- function(bin.size=1e3){
    if(!require(BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE)){
        stop("Please install BSgenome first by running:
> source(\"http://bioconductor.org/biocLite.R\")
> biocLite(\"BSgenome.Hsapiens.UCSC.hg19\")
")
    }
    seql.1.22 = seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[paste0("chr",1:22)]
    fragment.chr <- function(chr.i){
        starts = as.integer(seq(1,seql.1.22[chr.i], bin.size))
        ends = as.integer(starts + bin.size - 1)
        ends[length(ends)] = as.integer(seql.1.22[chr.i])
        data.frame(chr=chr.i,start = starts, end= ends)
    }
    plyr::ldply(lapply(1:22, fragment.chr), identity)
}
