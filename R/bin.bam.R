##' Count the number of reads from a BAM file in specified bins. By default the output bin counts are writtent into a file, which is then automatically compressed and indexed.
##'
##' By default, the function tries to check that the bin definition and the BAM file have compatible chromosome names. For example, 'chr1' in both or '1' in both. Moreover it will stop if no reads are found in the BAM file. To switch off these checks use 'no.checks=TRUE'.
##' @title Get read counts from BAM file
##' @param bam.file the BAM file or a list of BAM files (e.g. reverse and forward)
##' @param bin.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param outfile.prefix the prefix of the name of the output file. The suffix '.bgz' will
##' be appended to this name prefix after compression. In practice this is present in the output of 'init.filenames' function.
##' @param appendIndex.outfile if TRUE (default), the results will be appended regularly on
##' the output file which will be ultimately compressed and indexed. This is recommend when a large number
##' of bins are analyzed. If FALSE, a data.frame with the bin counts will be returned and no
##' file is created.
##' @param proper if TRUE (default), reads with properly mapped pairs are counted. If
##' FALSE, reads with improper mapping are counted.
##' @param map.quality the minimum mapping quality (Phred score) for the reads to count. Default
##' is 30.
##' @param chunk.size the number of bins to analyze at a time (for memory optimization).
##' Default is 10 000. Reduce this number if memory problems arise.
##' @param check.chr.name if TRUE (default), the function will try to check that the
##' definition of chromosome (e.g. '1' vs 'chr1') are consistent between the bin
##' definition and the BAM file. If FALSE, the analysis will continue either way.
##' @param no.checks if TRUE, won't stop if there are no reads or the chromosome names are inconsistent. Default is FALSE.
##' @param bai.file the index file of the BAM. If NULL it will be guessed.
##' @return a list:
##' \item{bc}{the final output file name if 'appendIndex.outfile' was TRUE; a data.frame with
##' the bin counts if not.}
##' \item{nb.reads}{the number of reads counted.}
##' @author Jean Monlong
##' @export
bin.bam <- function(bam.file, bin.df, outfile.prefix = NULL, appendIndex.outfile = TRUE,
                    proper = TRUE, map.quality = 30, chunk.size = 10000, check.chr.name = TRUE, no.checks=FALSE, bai.file=NULL) {
  if (is.null(outfile.prefix) & appendIndex.outfile) {
    stop("If 'appendIndex.outfile' is TRUE, please provide 'outfile.prefix'.")
  }
  if(!all(c("chr","start","end") %in% colnames(bin.df))){
    stop("Missing column in 'bin.df'. 'chr', 'start' and 'end' are required.")
  }
  if(any(!file.exists(bam.file))){
    stop("'bam.file' (",bam.file," file not found.")
  }

  bin.df = bin.df[order(as.character(bin.df$chr), bin.df$start),]
  bin.df$chunk = rep(1:ceiling(nrow(bin.df)/chunk.size), each = chunk.size)[1:nrow(bin.df)]
  bin.df$start = as.integer(bin.df$start)
  bin.df$end = as.integer(bin.df$end)

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
  
  if(length(bam.file)==1){
    binBam.single <- function(df) {
      gr.o = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start = start,
                                                                   end = end)))
      param = Rsamtools::ScanBamParam(which = gr.o, what = c("mapq"), flag = Rsamtools::scanBamFlag(isProperPair = proper,
                                                                                                    isDuplicate = FALSE, isNotPassingQualityControls = FALSE, isUnmappedQuery = FALSE))
      bam = Rsamtools::scanBam(bam.file, index = bai.file, param = param)
      unlist(lapply(bam, function(e) sum(unlist(e) > map.quality)))
    }
  } else {
    binBam.single <- function(df) {
      gr.o = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start = start,
                                                                   end = end)))
      param = Rsamtools::ScanBamParam(which = gr.o, what = c("mapq"), flag = Rsamtools::scanBamFlag(isProperPair = proper,
                                                                                                    isDuplicate = FALSE, isNotPassingQualityControls = FALSE, isUnmappedQuery = FALSE))
      bam.l = lapply(1:length(bam.file), function(ii){
        bam = Rsamtools::scanBam(bam.file[ii], index = bai.file[ii], param = param)
        unlist(lapply(bam, function(e) sum(unlist(e) > map.quality)))
      })
      do.call("+", bam.l)
    }
  }
  
  if (check.chr.name) {
    ## Check if headers
    bam.headers = Rsamtools::scanBamHeader(bam.file[1])
    if(length(bam.headers[[1]])>0){
      if(!any(paste0("SN:",unique(bin.df$chr)) %in% unlist(lapply(bam.headers[[1]][[2]],"[",1)))) {
        if (!grepl("chr", bin.df$chr[1])) {
          bin.df$chr = paste("chr", bin.df$chr, sep = "")
        } else {
          bin.df$chr = gsub("chr", "", bin.df$chr)
        }
      }
      if(!any(paste0("SN:",unique(bin.df$chr)) %in% unlist(lapply(bam.headers[[1]][[2]],"[",1))) & !no.checks) {
        stop("Couldn't guess if chr 1 is defined as '1' or 'chr1'.\nCheck manually and/or switch off option 'check.chr.name'.")
      }
    } else {
      ## Is it 'chr1' or '1': try with 20 random bins; if no read try with other
      bc.chrTest = binBam.single(bin.df[sample(1:nrow(bin.df), min(nrow(bin.df), 20)), ])
      if (all(bc.chrTest == 0)) {
        if (!grepl("chr", bin.df$chr[1])) {
          bin.df$chr = paste("chr", bin.df$chr, sep = "")
        } else {
          bin.df$chr = gsub("chr", "", bin.df$chr)
        }
        bc.chrTest = binBam.single(bin.df[sample(1:nrow(bin.df), min(nrow(bin.df), 20)), ])
        if (all(bc.chrTest == 0) & !no.checks) {
          stop("Couldn't guess if chr 1 is defined as '1' or 'chr1'.\nCheck manually and/or switch off option 'check.chr.name'.")
        }
      }
    }
  }

  binBam.chunk <- function(df) {
    ch.nb = as.numeric(df$chunk[1])
    df = df[, c("chr", "start", "end")]
    df$bc = binBam.single(df)
    if (!is.null(outfile.prefix)) {
      df$chunk = NULL
      utils::write.table(df, file = outfile.prefix, quote = FALSE, row.names = FALSE,
                  sep = "\t", append = ch.nb > 1, col.names = ch.nb == 1)
      return(data.frame(chunk = ch.nb, nb.reads = sum(as.numeric(df$bc), na.rm = TRUE)))
    } else {
      return(df)
    }
  }

  bc.df = lapply(unique(bin.df$chunk), function(chunk.i){
    binBam.chunk(bin.df[which(bin.df$chunk==chunk.i),])
  })
  bc.df = as.data.frame(data.table::rbindlist(bc.df))

  final.file = outfile.prefix
  if (appendIndex.outfile) {
    final.file = paste(final.file, ".bgz", sep = "")
    Rsamtools::bgzip(outfile.prefix, dest = final.file, overwrite = TRUE)
    file.remove(outfile.prefix)
    Rsamtools::indexTabix(final.file, format = "bed", skip=1L)
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
