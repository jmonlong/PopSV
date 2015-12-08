##' Correct GC bias by fitting a LOESS model. The bin counts are then corrected
##' with a continuous normalization factor derived from the predicted bin count
##' and global average bin count.
##' @title Correct GC bias
##' @param bc.f the name of the file with the bin count for a particular sample OR
##' a data.frame with 4 columns: 'chr', 'start', 'end' and 'bc'.
##' @param gc.df a data.frame with the bin definition and GC content information.
##' Column 'GCcontent' is required and can be obtained thanks to 'getGC.hp19' for
##' example. 
##' @param outfile.prefix the prefix for the output file name. The suffix '.bgz' will
##' be appended to this name prefix after compression. 
##' @param appendIndex.outfile if TRUE (default), the results will be saved on the
##' output file which will be ultimately compressed/indexed. This is recommend when a large number
##' of bins are analyzed. If FALSE, a data.frame with the corrected bin counts will
##' be returned and no file is created.
##' @return the name of the file with the corrected bin counts OR a new data.frame
##' with corrected bin counts.
##' @author Jean Monlong
##' @export
correct.GC <- function(bc.f, gc.df, outfile.prefix, appendIndex.outfile = TRUE) {
    gc.class = cut(gc.df$GCcontent, breaks = seq(0, 1, 0.02), include.lowest = TRUE)
    samp.ii = unlist(tapply(1:length(gc.class), gc.class, function(e) e[sample(1:length(e), 
        min(c(length(e), 500)))]))
    bc.df = read.table(bc.f, as.is = TRUE, header = TRUE)
    bc.df$GCcontent = gc.df$GCcontent
    lo = loess(bc ~ GCcontent, data = bc.df[samp.ii, ])
    bc.df$bc = mean(bc.df$bc, na.rm = TRUE) * bc.df$bc/predict(lo, newdata = bc.df)
    bc.df$bc = round(bc.df$bc, digits = 2)
    if (any(bc.df$bc < 0, na.rm = TRUE)) 
        bc.df$bc[bc.df$bc < 0] = 0
    bc.df$GCcontent = NULL
    
    if (appendIndex.outfile) {
        write.table(bc.df, file = outfile.prefix, quote = FALSE, row.names = FALSE, 
            sep = "\t")
        final.file = paste(outfile.prefix, ".bgz", sep = "")
        Rsamtools::bgzip(outfile.prefix, dest = final.file, overwrite = TRUE)
        file.remove(outfile.prefix)
        Rsamtools::indexTabix(final.file, format = "bed")
        return(final.file)
    } else {
        return(bc.df)
    }
    
} 
