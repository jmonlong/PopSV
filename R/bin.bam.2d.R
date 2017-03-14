##' Count the number of reads from a BAM file joining pairs of bins. By default the output bin counts are writtent into a file, which is then automatically compressed and indexed.
##'
##' By default, the function tries to check that the bin definition and the BAM file have compatible chromosome names. For example, 'chr1' in both or '1' in both. Moreover it will stop if no reads are found in the BAM file. To switch off these checks use 'no.checks=TRUE'.
##' @title Get read counts between pairs of bins from BAM file
##' @param bam.file the BAM file
##' @param bins.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param outfile.prefix the prefix of the name of the output file. The suffix '.bgz' will
##' be appended to this name prefix after compression. In practice this is present in the output of 'init.filenames' function.
##' @param appendIndex.outfile if TRUE (default), the results will be appended regularly on
##' the output file which will be ultimately compressed and indexed. This is recommend when a large number
##' of bins are analyzed. If FALSE, a data.frame with the bin counts will be returned and no
##' file is created.
##' @param chunk.size the number of bins to analyze at a time (for memory optimization).
##' Default is 1000. Reduce this number if memory problems arise.
##' @param check.chr.name if TRUE (default), the function will try to check that the
##' definition of chromosome (e.g. '1' vs 'chr1') are consistent between the bin
##' definition and the BAM file. If FALSE, the analysis will continue either way.
##' @param no.checks if TRUE, won't stop if there are no reads or the chromosome names are inconsistent. Default is FALSE.
##' @param bai.file the index file of the BAM. If NULL it will be guessed.
##' @return the final output file name if 'appendIndex.outfile' was TRUE; a data.frame with the bin counts if not.
##' @author Jean Monlong
##' @export
bin.bam.2d <- function(bam.file, bins.df, outfile.prefix = NULL, appendIndex.outfile = TRUE,
                    chunk.size = 1000, check.chr.name = TRUE, no.checks=FALSE, bai.file=NULL) {
  if (is.null(outfile.prefix) & appendIndex.outfile) {
    stop("If 'appendIndex.outfile' is TRUE, please provide 'outfile.prefix'.")
  }
  if(!all(c("chr","start","end") %in% colnames(bins.df))){
    stop("Missing column in 'bins.df'. 'chr', 'start' and 'end' are required.")
  }
  if(any(!file.exists(bam.file))){
    stop("'bam.file' (",bam.file," file not found.")
  }

  bins.df = bins.df[order(as.character(bins.df$chr), bins.df$start),]
  bins.df$chunk = as.numeric(cut(1:nrow(bins.df), nrow(bins.df)/chunk.size))
  bins.df$chrF = factor(bins.df$chr, levels=unique(bins.df$chr))

  if(is.null(bai.file)){
    bai.file = sub("bam$", "bai", bam.file, perl = TRUE)
  } else if(length(bai.file) != length(bam.file)){
    stop("Length of 'bam.file' different than length of 'bai.file'")
  }
  if (any(!file.exists(bai.file))) {
    bai.file = paste0(bam.file, ".bai")
    if (any(!file.exists(bai.file))) {
      stop("Index file is missing (neither '.bai' nor '.bam.bai').")
    }
  }
  
  count2D <- function(ch.i, bins.df, bam.file){
    bins.i = bins.df[which(bins.df$chunk==ch.i),]
    param = Rsamtools::ScanBamParam(which=GenomicRanges::makeGRangesFromDataFrame(bins.i),
                         what = c("mrnm","mpos"),
                         flag = Rsamtools::scanBamFlag(isPaired=TRUE,isProperPair=FALSE,isDuplicate=FALSE,isNotPassingQualityControls=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE)
                         )
    bam = Rsamtools::scanBam(bam.file,param=param)
    non.na.bin = which(unlist(lapply(bam, function(l)length(l[[1]])>0)))
    bam.df = lapply(non.na.bin, function(ii)data.frame(bin.i=ii, as.data.frame(bam[[ii]])))
    bam.df = as.data.frame(data.table::rbindlist(bam.df))
    bam.df$mrnmF = factor(as.character(bam.df$mrnm), levels=levels(bins.df$chrF))
    bam.df = cbind(bins.i[bam.df$bin.i,], bam.df)
    bam.df = bam.df[which(as.numeric(bam.df$mrnmF)>as.numeric(bam.df$chrF) |
                          (as.numeric(bam.df$mrnmF)>=as.numeric(bam.df$chrF) & bam.df$mpos>bam.df$end)),]
    if(nrow(bam.df)==0) return(NULL)
    bins.2.df = bins.df[which(bins.df$chunk>=ch.i),]
    bins.2 = GenomicRanges::makeGRangesFromDataFrame(bins.2.df)
    bam.2 = GenomicRanges::GRanges(bam.df$mrnm, IRanges::IRanges(bam.df$mpos, width=1))
    ol2 = GenomicRanges::findOverlaps(bam.2, bins.2)
    bam.df$start2 = bam.df$chr2 = NA
    bam.df$chr2[S4Vectors::queryHits(ol2)] = bins.2.df$chr[S4Vectors::subjectHits(ol2)]
    bam.df$start2[S4Vectors::queryHits(ol2)] = bins.2.df$start[S4Vectors::subjectHits(ol2)]
    bam.df = bam.df[,c("chr","start","end","chr2","start2")]
    bam.df$read = 1
    bam.df = stats::aggregate(read~chr+start+end+chr2+start2, data=bam.df, sum)
    bam.df = bam.df[order(as.character(bam.df$chr), bam.df$start),]    
    if (!is.null(outfile.prefix)) {
      utils::write.table(bam.df, file = outfile.prefix, quote = FALSE, row.names = FALSE,
                         sep = "\t", append = ch.i > 1, col.names = ch.i == 1)
      return(data.frame(chunk = ch.i, nb.reads = sum(as.numeric(bam.df$read), na.rm = TRUE)))
    } else {
      return(bam.df)
    }    
  }
  
  if (check.chr.name) {
    ## Check if headers
    bam.headers = Rsamtools::scanBamHeader(bam.file[1])
    if(length(bam.headers[[1]])>0){
      if(!any(paste0("SN:",unique(bins.df$chr)) %in% unlist(lapply(bam.headers[[1]][[2]],"[",1)))) {
        if (!grepl("chr", bins.df$chr[1])) {
          bins.df$chr = paste("chr", bins.df$chr, sep = "")
        } else {
          bins.df$chr = gsub("chr", "", bins.df$chr)
        }
      }
      if(!any(paste0("SN:",unique(bins.df$chr)) %in% unlist(lapply(bam.headers[[1]][[2]],"[",1))) & !no.checks) {
        stop("Couldn't guess if chr 1 is defined as '1' or 'chr1'.\nCheck manually and/or switch off option 'check.chr.name'.")
      }
    } 
  }

  bc.df = lapply(unique(bins.df$chunk), count2D, bins.df=bins.df, bam.file=bam.file)
  bc.df = as.data.frame(data.table::rbindlist(bc.df))

  final.file = outfile.prefix
  if (appendIndex.outfile) {
    final.file = paste(final.file, ".bgz", sep = "")
    Rsamtools::bgzip(outfile.prefix, dest = final.file, overwrite = TRUE)
    file.remove(outfile.prefix)
    Rsamtools::indexTabix(final.file, format = "bed")
  }
  if(is.null(outfile.prefix)) {
    bc.df$chunk = NULL
    res.l = list(bc = bc.df, nb.reads = sum(as.numeric(bc.df$bc), na.rm = TRUE))
  } else {
    res.l = list(bc = final.file, nb.reads = sum(as.numeric(bc.df$nb.reads), na.rm = TRUE))
  }
  if(res.l$nb.reads==0 & !no.checks){
    stop("no reads found with these parameters (binning, mapping) in:",bam.file,". Suspicious, maybe investigate.")
  }
  return(res.l)
}
