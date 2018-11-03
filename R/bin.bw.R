##' Count the number of reads from a BigWig file in specified bins. By default the output bin counts are writtent into a file, which is then automatically compressed and indexed.
##'
##' By default, the function tries to check that the bin definition and the BAM file have compatible chromosome names. For example, 'chr1' in both or '1' in both. Moreover it will stop if no reads are found in the BAM file. To switch off these checks use 'no.checks=TRUE'.
##' @title Get read counts from BigWig coverage file
##' @param bw.file the BigWig file
##' @param bin.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param outfile.prefix the prefix of the name of the output file. The suffix '.bgz' will
##' be appended to this name prefix after compression. In practice this is present in the output of 'init.filenames' function.
##' @param appendIndex.outfile if TRUE (default), the results will be appended regularly on
##' the output file which will be ultimately compressed and indexed. This is recommend when a large number
##' of bins are analyzed. If FALSE, a data.frame with the bin counts will be returned and no
##' file is created.
##' @param chunk.size the number of bins to analyze at a time (for memory optimization).
##' Default is 10 000. Reduce this number if memory problems arise.
##' @param check.chr.name if TRUE (default), the function will try to check that the
##' definition of chromosome (e.g. '1' vs 'chr1') are consistent between the bin
##' definition and the BAM file. If FALSE, the analysis will continue either way.
##' @param no.checks if TRUE, won't stop if there are no reads or the chromosome names are inconsistent. Default is FALSE.
##' @param read.length the read length, to normalize the counts. Default is 100.
##' @param fromSummaries should the count be derived from the BigWig summaries (much faster). Default is TRUE.
##' @return a list:
##' \item{bc}{the final output file name if 'appendIndex.outfile' was TRUE; a data.frame with
##' the bin counts if not.}
##' \item{nb.reads}{the number of reads counted.}
##' @author Jean Monlong
##' @export
bin.bw <- function(bw.file, bin.df, outfile.prefix = NULL, appendIndex.outfile = TRUE,
                    chunk.size = 10000, check.chr.name = TRUE, no.checks=FALSE, read.length=100, fromSummaries=TRUE) {
  if (is.null(outfile.prefix) & appendIndex.outfile) {
    stop("If 'appendIndex.outfile' is TRUE, please provide 'outfile.prefix'.")
  }
  if(!all(c("chr","start","end") %in% colnames(bin.df))){
    stop("Missing column in 'bin.df'. 'chr', 'start' and 'end' are required.")
  }
  if(any(!file.exists(bw.file))){
    stop("'bw.file' (",bw.file," file not found.")
  }

  ## Uglily appeases R checks
  subjectHits = NULL

  if(fromSummaries){
    binBam.single <- function(df) {
      gr.o = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start = start, end = end)))
      scores = summary(bwf, gr.o, defaultValue=0, type='mean')
      df$bc = as.data.frame(scores)$score * width(gr.o) / read.length
      df
    }
  } else {
    binBam.single <- function(df) {
      gr.o = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start = start, end = end)))
      scores = rtracklayer::import(bw.file, which=gr.o, format="BigWig", as="NumericList")
      df = as.data.frame(S4Vectors::metadata(scores)$ranges)
      df = df[, 1:3]
      colnames(df)[1] = 'chr'
      df$chr = as.character(df$chr)
      df$bc = unlist(lapply(scores, sum)) / read.length
      df
    }
  }

  bwf <- rtracklayer::BigWigFile(bw.file)
  chr.levels = seqnames(seqinfo(bwf))

  if (check.chr.name) {
    if(all(grepl('chr', chr.levels))){
      if (!grepl("chr", bin.df$chr[1])) {
        bin.df$chr = paste("chr", bin.df$chr, sep = "")
      }
    } else {
      if (grepl("chr", bin.df$chr[1])) {
        bin.df$chr = gsub("chr", "", bin.df$chr)
      }
    }
  }

  bin.df = bin.df[order(as.character(bin.df$chr), bin.df$start),]
  bin.df$chunk = rep(1:ceiling(nrow(bin.df)/chunk.size), each = chunk.size)[1:nrow(bin.df)]

  binBam.chunk <- function(df) {
    ch.nb = as.numeric(df$chunk[1])
    df = df[, c("chr", "start", "end")]
    df = binBam.single(df)
    df = df[order(df$chr, df$start),]
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
    Rsamtools::indexTabix(final.file, format = "bed")
  }
  if(is.null(outfile.prefix)) {
    bc.df$chunk = NULL
    res.l = list(bc = bc.df, nb.reads = sum(as.numeric(bc.df$bc), na.rm = TRUE))
  } else {
    res.l = list(bc = final.file, nb.reads = sum(as.numeric(bc.df$nb.reads), na.rm = TRUE))
  }
  if(res.l$nb.reads==0 & !no.checks){
    stop("no reads found with these parameters (binning, mapping) in:",bw.file,". Suspicious, maybe investigate.")
  }
  return(res.l)
}
