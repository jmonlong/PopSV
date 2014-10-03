##' Using \code{Rsamtools} functions to manipulate tabix file, the indexed BED file is read
##' and only regions defined by the desired subset are retrieved. 
##' @title Retrieve a subset of an indexed BED file. 
##' @param file the name of the file which the data are to be read from.
##' Should be compressed with |code{bgzip} and indexed by tabix algorithm.
##' @param subset a data.frame or GRanges object with the regions to subset from.
##' If a data.frame, it must have columns named 'chr', 'start' and 'end'.
##' @param col.names The header for the output object. Optional.
##' @param as.is controls the conversion of columns. Default is TRUE, i.e. columns are
##' converted into 'character', 'numeric', 'integer', etc. If FALSE, some columns might
##' be converted into factors (e.g. 'character' columns).
##' @return a data.frame with the retrieved BED information.
##' @author Jean Monlong
read.bedix <- function(file,subset, col.names=NULL, as.is=TRUE){

    if(is.data.frame(subset)){
        subset = with(subset, GenomicRanges::GRanges(chr, IRanges(start, end)))
    } else if(class(subset)!="GRanges"){
        stop("'subset' must be a data.frame or a GRanges object.")
    }

    bed = unlist(Rsamtools::scanTabix(file,param=subset))
    ncol = length(strsplit(bed[1],"\t")[[1]])
    bed.df = matrix(unlist(strsplit(bed,"\t")), length(bed), ncol, byrow=TRUE)
    bed.df = as.data.frame(bed.df, stringsAsFactors=FALSE)
    if(!is.null(col.names)){
        colnames(bed.df) = col.names
    } else {
      colnames(bed.df) = as.character(read.table(file, nrows=1, as.is=TRUE))
    }
    for(ii in 1:ncol(bed.df)) bed.df[,ii] = type.convert(bed.df[,ii], as.is=as.is)
    
    return(bed.df)
}
