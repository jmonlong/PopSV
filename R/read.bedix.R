##' Using \code{Rsamtools} functions to manipulate tabix file, the indexed BED file is read
##' and only regions defined by the desired subset are retrieved. 
##' @title Retrieve a subset of an indexed BED file. 
##' @param file the name of the file which the data are to be read from.
##' Should be compressed with |code{bgzip} and indexed by tabix algorithm.
##' @param subset.reg a data.frame or GRanges object with the regions to subset from.
##' If a data.frame, it must have columns named 'chr', 'start' and 'end'.
##' @param header Should the first line be used as headers. Default is TRUE.
##' @param as.is controls the conversion of columns. Default is TRUE, i.e. columns are
##' converted into 'character', 'numeric', 'integer', etc. If FALSE, some columns might
##' be converted into factors (e.g. 'character' columns).
##' @return a data.frame with the retrieved BED information.
##' @author Jean Monlong
##' @export
##' @import data.table
read.bedix <- function(file, subset.reg=NULL, header=TRUE, as.is = TRUE) {

  if(!file.exists(file)){
    file = paste0(file, ".bgz")
  }
  if(!file.exists(file)){
    stop(file, "Input file not found (with and without .bgz extension).")
  }
  if(!file.exists(paste0(file,".tbi"))){
    stop(paste0(file,".tbi"),"Index file not found.")
  }
  if(is.null(subset.reg)){
    return(read.table(file, as.is=as.is, header=header))
  }
  if (is.data.frame(subset.reg)) {
    if(!all(c("chr","start","end") %in% colnames(subset.reg))){
      stop("Missing column in 'subset.reg'. 'chr', 'start' and 'end' are required.")
    }
    subset.reg = with(subset.reg, GenomicRanges::GRanges(chr, IRanges::IRanges(start,end)))
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
    if (header) {
      data.table::setnames(bed, as.character(read.table(file, nrows = 1, as.is = TRUE)))
    }
    bed = bed[, lapply(.SD, function(ee)type.convert(as.character(ee), as.is=TRUE))]
    bed = as.data.frame(bed)
    return(bed)
  }
  
  if (length(subset.reg) > 10000) {
    chunks = cut(1:length(subset.reg), ceiling(length(subset.reg)/10000))
    bed.df = plyr::ldply(levels(chunks), function(ch.id){
      read.chunk(subset.reg[which(chunks == ch.id)])
    })
  } else {
    bed.df = read.chunk(subset.reg)
  }
  bed.df = bed.df[order(bed.df$chr, bed.df$start),]
  return(bed.df)
} 
