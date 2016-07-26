##' Count the number of discordant reads from a BAM file in specified bins. By default the output bin counts are writtent into a file, which is then automatically compressed and indexed.
##'
##' By default, the function tries to check that the bin definition and the BAM file have compatible chromosome names. For example, 'chr1' in both or '1' in both. Moreover it will stop if no reads are found in the BAM file. To switch off these checks use 'no.checks=TRUE'.
##' @title Get discordant read counts from BAM file
##' @param bam.file the BAM file
##' @param bin.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param outfile.prefix the prefix of the name of the output file. The suffix '.bgz' will
##' be appended to this name prefix after compression. In practice this is present in the output of 'init.filenames' function.
##' @param appendIndex.outfile if TRUE (default), the results will be appended regularly on
##' the output file which will be ultimately compressed and indexed. This is recommend when a large number
##' of bins are analyzed. If FALSE, a data.frame with the bin counts will be returned and no
##' file is created.
##' @param min.weight the minimum weight (when a discordant read is isolated). Default is 0.1.
##' @param weight.dist the distance from its neighbor where a read is considered isolated in the weighting process. Default is 50 (bp).
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
bin.bam.disc <- function(bam.file, bin.df, outfile.prefix = NULL, appendIndex.outfile = TRUE,
                           min.weight=.1, weight.dist = 50, chunk.size = 10000, check.chr.name = TRUE, no.checks=FALSE, bai.file=NULL) {
  if (is.null(outfile.prefix) & appendIndex.outfile) {
    stop("If 'appendIndex.outfile' is TRUE, please provide 'outfile.prefix'.")
  }
  if(!all(c("chr","start","end") %in% colnames(bin.df))){
    stop("Missing column in 'bin.df'. 'chr', 'start' and 'end' are required.")
  }
  if(!file.exists(bam.file)){
    stop("'bam.file' (",bam.file," file not found.")
  }

  bin.df = bin.df[order(as.character(bin.df$chr), bin.df$start),]
  bin.df$chunk = rep(1:ceiling(nrow(bin.df)/chunk.size), each = chunk.size)[1:nrow(bin.df)]

  if(is.null(bai.file)){
    bai.file = sub("bam$", "bai", bam.file, perl = TRUE)
  }
  if (!file.exists(bai.file)) {
    bai.file = paste0(bam.file, ".bai")
    if (!file.exists(bai.file)) {
      stop("Index file is missing (neither '.bai' nor '.bam.bai').")
    }
  }

  flagCigar <- function(cigar){
    cigar = ifelse(is.na(cigar),"NAS",cigar)
    cigar.l = gsub("[0-9]+M","0S",as.character(cigar))
    cigar.l = strsplit(cigar.l, "[SIDHP]",perl=TRUE)
    cigar.l = suppressWarnings(lapply(cigar.l, as.numeric))
    sapply(cigar.l,max)
  }
  sameOrientation <- function(flag){
    ##proper = floor(flag / 2) %% 2 == 1
    rrev = floor(flag / 2**4) %% 2
    mrev = floor(flag / 2**5) %% 2
    rrev == mrev
  }
  Rcpp::cppFunction('
NumericVector DPLOOP(NumericVector p, NumericVector id)
{
  int n = p.size() / 2;
  NumericVector res(n);
  int maxp = max(p);
  for(int i = 0; i < res.size(); ++i)
    {
      res[i] = maxp;
    }
  for(int i = 0; i < p.size(); ++i)
    {
      int is = i-2;
      if(is < 0){
	is = 0;
      }
      int ie = i+2;
      if(ie > p.size()){
	ie = p.size();
      }
      for(int j = is; j < ie; ++j)
	{
	  if(id[j] != id[i] & res[id[i]] > abs(p[j]-p[i]))
	    {
	      res[id[i]] = abs(p[j]-p[i]);
	    }
	}
    }
  return wrap(res);
}
')
  distpairC <- function(x, y){
    dfv = data.frame(p=c(x,x), id=0:(length(x)-1))
    dfv = dfv[order(dfv$p),]
    DPLOOP(dfv$p, dfv$id)
  }
  binBam.single <- function(bins) {
    gr.o = with(bins, GenomicRanges::GRanges(chr, IRanges::IRanges(start = start, end = end)))
    param = Rsamtools::ScanBamParam(which = gr.o, what = c("rname","pos","flag","isize","mrnm","mpos","cigar","mapq"))
    reads.l = Rsamtools::scanBam(bam.file, index = bai.file, param = param)
    insert.size = unlist(lapply(reads.l, function(e)abs(e$isize)))
    is.med = median(insert.size, na.rm=TRUE)
    is.dev = abs(insert.size-is.med)
    is.q = quantile(is.dev, na.rm=TRUE, probs=.99)
    scores.l = lapply(reads.l, function(ll){
      if(length(ll[[1]])==0){
        return(0)
      }
      is.ol = abs(abs(ll$isize)-is.med) > is.q
      dbp=flagCigar(ll$cigar)
      same.or=sameOrientation(ll$flag)
      disc=is.ol | same.or | is.na(ll$isize) | dbp > 50 ## Discordant mapping definition (recipe)
      disc = which(disc & !is.na(ll$pos))
      if(length(disc)==0){
        return(0)
      }
      dw = distpairC(ll$pos[disc], ll$mpos[disc])
      dw = 1-(1-min.weight)*as.numeric(dw)/weight.dist
      dw = ifelse(dw<min.weight, min.weight, dw)
      return(sum(dw))
    })
    return(unlist(scores.l))
  }

  if (check.chr.name) {
    ## Check if headers
    bam.headers = Rsamtools::scanBamHeader(bam.file)
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
      write.table(df, file = outfile.prefix, quote = FALSE, row.names = FALSE,
                  sep = "\t", append = ch.nb > 1, col.names = ch.nb == 1)
      return(data.frame(chunk = ch.nb, nb.reads = sum(as.numeric(df$bc), na.rm = TRUE)))
    } else {
      return(df)
    }
  }

  bc.df = lapply(unique(bin.df$chunk), function(chunk.i){
    binBam.chunk(bin.df[which(bin.df$chunk==chunk.i),])
  })
  bc.df = do.call(rbind, bc.df)

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
