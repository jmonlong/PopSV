##' Correct GC bias by fitting a LOESS model. The bin counts are then corrected
##' with a continuous normalization factor derived from the predicted bin count
##' and global average bin count.
##' @title Correct GC bias
##' @param bc.f the name of the file with the bin count for a particular sample OR
##' a data.frame with 4 columns: 'chr', 'start', 'end' and 'bc'.
##' @param gc.df a data.frame with the bin definition and GC content information.
##' Column 'GCcontent' is required and can be obtained thanks to 'getGC.hg19' for
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
    if(!is.data.frame(bc.f) & length(bc.f)==1 & is.character(bc.f)){
      bc.df = utils::read.table(bc.f, as.is = TRUE, header = TRUE)
    } else {
      bc.df = bc.f
    }
    gc.df$start = as.integer(gc.df$start)
    gc.df$end = as.integer(gc.df$end)
    bc.df$start = as.integer(bc.df$start)
    bc.df$end = as.integer(bc.df$end)
    if(any(bc.df$chr != gc.df$chr) | any(bc.df$start != gc.df$start)){
      bc.df = merge(bc.df, gc.df[,c('chr','start','end','GCcontent')])
      bc.df = bc.df[order(as.character(bc.df$chr), bc.df$start),]
    } else {
      bc.df$GCcontent = gc.df$GCcontent
    }
    lo = stats::loess(bc ~ GCcontent, data = bc.df[samp.ii, ])
    bc.df$bc = mean(bc.df$bc, na.rm = TRUE) * bc.df$bc/stats::predict(lo, newdata = bc.df)
    bc.df$bc = round(bc.df$bc, digits = 2)
    if (any(bc.df$bc < 0, na.rm = TRUE))
        bc.df$bc[bc.df$bc < 0] = 0
    bc.df$GCcontent = NULL

    if (appendIndex.outfile) {
        utils::write.table(bc.df, file = outfile.prefix, quote = FALSE, row.names = FALSE,
            sep = "\t")
        final.file = paste(outfile.prefix, ".bgz", sep = "")
        Rsamtools::bgzip(outfile.prefix, dest = final.file, overwrite = TRUE)
        file.remove(outfile.prefix)
        Rsamtools::indexTabix(final.file, format = "bed", skip=1L)
        return(final.file)
    } else {
        return(bc.df)
    }

}
