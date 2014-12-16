##' Compress and index BED-like files. The files need to be ordered by genomic position and the first three columns must define chromosome location, start and end position.
##' @title Compress and index BED-like files
##' @param files the names of the files to compress and index.
##' @param outprefix the prefix to use to name the output files. By default, the original file names.
##' @param rm.input Should the original input file be removed after compression/indexing ? Default is TRUE.
##' @param overwrite.out Should a compressed file be overwritten if it already exists ? Default is TRUE.
##' @return a character vector confirming the creation of the new files.
##' @author Jean Monlong
##' @export
comp.index.files <- function(files, outprefix=files, rm.input=TRUE, overwrite.out=TRUE){
    sapply(1:length(files), function(file.ii){
        final.file = paste(outprefix[file.ii],".bgz",sep="")
        Rsamtools::bgzip(files[file.ii], dest=final.file, overwrite=TRUE)
        if(rm.input) file.remove(files[file.ii])
        Rsamtools::indexTabix(final.file, format="bed")
        return(paste(final.file,"created and indexed."))
    })
}
