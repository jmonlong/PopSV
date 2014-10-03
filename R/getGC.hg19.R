##' Add a column with GC content to a data.frame with 'chr', 'start' and 'end' column.
##' Requires BSgenome for Hsapiens hg19; if not installed instructions to install it will
##' be displayed.
##' @title GC content computation for specific bins
##' @param bins.df a data.frame with bin information. Columns 'chr', 'start' and 'end'
##' are required.
##' @return the same input data.frame with an extra column 'GCcontent'.
##' @author Jean Monlong
##' @export
getGC.hg19 <- function(bins.df){
    if(!require(BSgenome.Hsapiens.UCSC.hg19)){
        stop("Please install BSgenome first by running:
> source(\"http://bioconductor.org/biocLite.R\")
> biocLite(\"BSgenome.Hsapiens.UCSC.hg19\")
")
    }
    bins.df$chunk = rep(1:ceiling(nrow(bins.df)/1e3),each=1e3)[1:nrow(bins.df)]
    addGC <- function(df){
        if(!grepl("chr",df$chr[1])){
            chrs = paste("chr",df$chr,sep="")
        } else {
            chrs = df$chr
        }
        seq.l = Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chrs, df$start, df$end)
        lf = Biostrings::letterFrequency(seq.l,letters=c("G","C"))
        df$GCcontent = rowSums(lf) / (df$end-df$start)
        subset(df, !is.na(GCcontent))
    }
    bins.df = dplyr::do(dplyr::group_by(bins.df,chunk),addGC(.))
    bins.df$chunk = NULL
    return(bins.df)
}
