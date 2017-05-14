##' Compress and index BED-like files. The first three columns must define chromosome location, start and end position.
##'
##' The files should be ordered by genomic position. If not 'reorder=TRUE' will force the files to be imported in R, reordered and written back.
##' @title Compress and index BED-like files
##' @param files the names of the files to compress and index.
##' @param outprefix the prefix to use when naming the output files. By default, the original file names.
##' @param rm.input should the original input file be removed after compression/indexing ? Default is TRUE.
##' @param overwrite.out should a compressed file be overwritten if it already exists ? Default is TRUE.
##' @param reorder should the files be read and reordered before compression/indexation. Default is FALSE.
##' @return a character vector confirming the creation of the new files.
##' @author Jean Monlong
##' @import data.table
##' @export
comp.index.files <- function(files, outprefix = files, rm.input = TRUE, overwrite.out = TRUE, reorder=FALSE) {
  chr = start = NULL ## Uglily appeases R checks
  if(any(!file.exists(files))){
    stop(files[which(!file.exists(files))], ": file not found")
  }
  sapply(1:length(files), function(file.ii) {
           if(reorder){
             dt = data.table::fread(files[file.ii])
             dt[, chr:= as.character(chr)]
             if(all(c("chr2", "start2") %in% colnames(dt))){
               data.table::setkey(dt, chr, start, chr2, start2)
             } else {
               data.table::setkey(dt, chr, start)
             }
             utils::write.table(dt, file=files[file.ii], quote=FALSE, row.names=FALSE, sep="\t")
           }
           final.file = paste(outprefix[file.ii], ".bgz", sep = "")
           Rsamtools::bgzip(files[file.ii], dest = final.file, overwrite = TRUE)
           if (rm.input)
             file.remove(files[file.ii])
           Rsamtools::indexTabix(final.file, format = "bed")
           ##message(paste(final.file, "created and indexed."))
           return(final.file)
         })
}
