##' Using \code{Rsamtools} functions to manipulate tabix file, the indexed BED file is read
##' and only regions defined by the desired subset are retrieved. 
##' @title Retrieve a subset of an indexed BED file. 
##' @param file the name of the file which the data are to be read from.
##' Should be compressed with |code{bgzip} and indexed by tabix algorithm.
##' @param subset.reg a data.frame or GRanges object with the regions to subset from.
##' If a data.frame, it must have columns named 'chr', 'start' and 'end'.
##' @param col.names The header for the output object. Optional.
##' @param as.is controls the conversion of columns. Default is TRUE, i.e. columns are
##' converted into 'character', 'numeric', 'integer', etc. If FALSE, some columns might
##' be converted into factors (e.g. 'character' columns).
##' @return a data.frame with the retrieved BED information.
##' @author Jean Monlong
##' @export
read.bedix <- function(file, subset.reg, col.names = NULL, as.is = TRUE) {
    
    if (is.data.frame(subset.reg)) {
        subset.reg = with(subset.reg, GenomicRanges::GRanges(chr, IRanges::IRanges(start, 
            end)))
    } else if (class(subset.reg) != "GRanges") {
        stop("'subset.reg' must be a data.frame or a GRanges object.")
    }

    subset.reg = subset.reg[order(as.character(GenomicRanges::seqnames(subset.reg)), GenomicRanges::start(subset.reg))]

    read.chunk <- function(gr){
      bed = tryCatch(unlist(Rsamtools::scanTabix(file, param = GenomicRanges::reduce(gr))), error = function(e) c())
      if (length(bed) == 0) {
        return(NULL)
      }
      ncol = length(strsplit(bed[1], "\t")[[1]])
      bed = matrix(unlist(strsplit(bed, "\t")), length(bed), ncol, byrow = TRUE)
      bed = data.table::data.table(bed)
      bed = bed[, lapply(.SD, type.convert, as.is=TRUE)]
      bed = as.data.frame(bed)
      ##bed = as.data.frame(bed, stringsAsFactors = FALSE)
      if (!is.null(col.names)) {
        colnames(bed) = col.names
      } else {
        colnames(bed) = as.character(read.table(file, nrows = 1, as.is = TRUE))
      }
      return(bed)
    }
    if (length(subset.reg) > 10000) {
      chunks = cut(1:length(subset.reg), ceiling(length(subset.reg)/10000))
      bed.df = plyr::ldply(levels(chunks), function(ch.id){
        cat(ch.id,"\n")
        read.chunk(subset.reg[which(chunks == ch.id)])
      })
    } else {
      bed.df = read.chunk(subset.reg)
    }
    
    return(bed.df)
} 
