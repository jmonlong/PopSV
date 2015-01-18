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
read.bedix <- function(file,subset.reg, col.names=NULL, as.is=TRUE){

    if(is.data.frame(subset.reg)){
        subset.reg = with(subset.reg, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    } else if(class(subset.reg)!="GRanges"){
        stop("'subset.reg' must be a data.frame or a GRanges object.")
    }

    subset.reg = GenomicRanges::reduce(subset.reg)
    
    bed = tryCatch(unlist(Rsamtools::scanTabix(file,param=subset.reg)),
        error=function(e)c())
    if(length(bed)==0){
        return(NULL)
    }
    gc() ## Not sure if needed
    ncol = length(strsplit(bed[1],"\t")[[1]])
    if(length(bed)>1e4){
      bed.df = matrix(NA, length(bed), ncol)
      chunks = cut(1:length(bed), ceiling(length(bed)/1e4))
      for(ch.id in levels(chunks)){
        ch.ii = which(chunks==ch.id)
        bed.df[ch.ii,] = matrix(unlist(strsplit(bed[ch.ii],"\t")), length(ch.ii), ncol, byrow=TRUE)
      }
    } else {
      bed.df = matrix(unlist(strsplit(bed,"\t")), length(bed), ncol, byrow=TRUE)
    }
    rm(bed)
    bed.df = as.data.frame(bed.df, stringsAsFactors=FALSE)
    if(!is.null(col.names)){
        colnames(bed.df) = col.names
    } else {
      colnames(bed.df) = as.character(read.table(file, nrows=1, as.is=TRUE))
    }
    gc()
    col.classes = c("character",rep("integer",2), rep("numeric",ncol-3))
    for(ii in 1:ncol(bed.df)){
      class(bed.df[,ii]) = col.classes[ii]
      gc()
    }
    
    return(bed.df)
}
