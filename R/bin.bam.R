##' Count the number of reads from a BAM file in specified bins. 
##' @title Get read counts from BAM file
##' @param bam.file the BAM file
##' @param bin.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param outfile.prefix the prefix of the name of the output file. The suffix '.bgz' will
##' be appended to this name prefix after compression. 
##' @param appendIndex.outfile if TRUE (default), the results will be appended regularly on
##' the output file which will be ultimately indexed. This is recommend when a large number
##' of bins are analyzed. If FALSE, a data.frame with the bin counts will be returned and no
##' file are created.
##' @param proper if TRUE (default), reads with properly mapped pairs are counted. If
##' FALSE, reads with improper mapping are counted.
##' @param map.quality the minimum mapping quality (PHRED) for the reads to count. Default
##' is 30.
##' @param chunk.size the number of bins to analyze at a time (for memory optimization).
##' Default is 10 000. Reduce this number if memory problems arise.
##' @param check.chr.name if TRUE (default), the function will try to check that the
##' definition of chromosome (e.g. '1' vs 'chr1') are consistent between the bin
##' definition and the BAM file. If FALSE, the analysis will continue either way.
##' @return a list:
##' \item{bc}{the final output file name if 'appendIndex.outfile' was TRUE; a data.frame with
##' the bin counts if not.}
##' \item{nb.reads}{the number of reads counted.}
##' @author Jean Monlong
##' @export
bin.bam <- function(bam.file,bin.df,outfile.prefix=NULL, appendIndex.outfile=TRUE,proper=TRUE,map.quality=30, chunk.size=1e4, check.chr.name=TRUE){
    if(is.null(outfile.prefix) & appendIndex.outfile){
        stop("If 'appendIndex.outfile' is TRUE, please provide 'outfile.prefix'.")
    }
    
    bin.df = dplyr::arrange(bin.df, chr, start)
    bin.df$chunk = rep(1:ceiling(nrow(bin.df)/chunk.size),each=chunk.size)[1:nrow(bin.df)]

    bai.file = sub("bam$","bai",bam.file,perl=TRUE)
    if(!file.exists(bai.file)){
        bai.file = paste0(bam.file,".bai")
        if(!file.exists(bai.file)){
            stop("Index file is missing (neither '.bai' nor '.bam.bai').")
        }
    }

    binBam.single <- function(df){
        gr.o = with(df,GenomicRanges::GRanges(chr,IRanges::IRanges(start=start,end=end)))
        param = Rsamtools::ScanBamParam(which=gr.o,
            what = c("mapq"),
            flag = Rsamtools::scanBamFlag(isProperPair=proper,isDuplicate=FALSE,isNotPassingQualityControls=FALSE,isUnmappedQuery=FALSE)
                                        )
        bam = Rsamtools::scanBam(bam.file, index=bai.file,param=param)
        unlist(lapply(bam,function(e)sum(unlist(e)>map.quality)))
    }

    if(check.chr.name){
        ## Is it "chr1" or "1": try with 10 random bins; if no read try with other
        bc.chrTest = binBam.single(bin.df[sample(1:nrow(bin.df),min(nrow(bin.df),10)),])
        if(all(bc.chrTest==0)){
            if(!grepl("chr",bin.df$chr[1])){
                bin.df$chr = paste("chr",bin.df$chr,sep="")
            } else {
                bin.df$chr = gsub("chr","",bin.df$chr)
            }
            bc.chrTest = binBam.single(bin.df[sample(1:nrow(bin.df),min(nrow(bin.df),10)),])
            if(all(bc.chrTest==0)){
                stop("Couldn't guess if chr 1 is defined as '1' or 'chr1'.
Check manually and/or switch off option 'check.chr.name'.")
            }
        }
    }

    binBam.chunk <- function(df){
        ch.nb = as.numeric(df$chunk[1])
        df = df[,c("chr","start","end")]
        df$bc = binBam.single(df)
        if(appendIndex.outfile & !is.null(outfile.prefix)){
            df$chunk = NULL
            write.table(df, file=outfile.prefix, quote=FALSE, row.names=FALSE, sep="\t", append=ch.nb>1, col.names=ch.nb==1)
            return(data.frame(chunk=ch.nb, nb.reads=sum(df$bc, na.rm=TRUE)))
        } else {
            return(df)
        }
    }

    bc.df = dplyr::do(dplyr::group_by(bin.df,chunk),binBam.chunk(.))

    if(appendIndex.outfile & !is.null(outfile.prefix)){
        final.file = paste(outfile.prefix,".bgz",sep="")
        Rsamtools::bgzip(outfile.prefix, dest=final.file, overwrite=TRUE)
        file.remove(outfile.prefix)
        Rsamtools::indexTabix(final.file, format="bed")
        return(list(bc=final.file, nb.reads=sum(bc.df$nb.reads, na.rm=TRUE)))
    } else {
        bc.df$chunk = NULL
        return(list(bc=bc.df, nb.reads=sum(bc.df$bc, na.rm=TRUE)))
    }    
}
