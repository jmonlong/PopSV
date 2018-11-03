##' Add a column with GC content to a data.frame with 'chr', 'start' and 'end' column.
##' Requires BSgenome for Hsapiens hg19; if not instructions to install it will
##' be displayed.
##' @title GC content computation for specific bins
##' @param bins.df a data.frame with bin information. Columns 'chr', 'start' and 'end'
##' are required.
##' @param genome the BSgenome object with the genome. Usually BSgenome.Hsapiens.UCSC.hg19 or BSgenome.Hsapiens.UCSC.hg38.
##' @return the same input data.frame with an extra column 'GCcontent'.
##' @author Jean Monlong
##' @export
getGC <- function(bins.df, genome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {
    if(!all(c("chr","start","end") %in% colnames(bins.df))){
    stop("Missing column in 'bin.df'. 'chr', 'start' and 'end' are required.")
  }
    bins.df$chunk = rep(1:ceiling(nrow(bins.df)/1000), each = 1000)[1:nrow(bins.df)]
    addGC <- function(df) {
        if (!grepl("chr", df$chr[1])) {
            chrs = paste("chr", df$chr, sep = "")
        } else {
            chrs = df$chr
        }
        seq.l = Biostrings::getSeq(genome, chrs, df$start, df$end)
        lf = Biostrings::letterFrequency(seq.l, letters = c("G", "C"))
        df$GCcontent = rowSums(lf)/(df$end - df$start + 1)
        df[which(!is.na(df$GCcontent)), ]
    }
    ## bins.df = dplyr::do(dplyr::group_by(bins.df,chunk),addGC(.))
    chunk = . = NULL  ## Uglily appease R checks
    bins.df = dplyr::do(dplyr::group_by(bins.df, chunk), addGC(.))
    bins.df$chunk = NULL
    return(as.data.frame(bins.df))
}
