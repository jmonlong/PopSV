##' Compare the pair mapping between regions of the genome. The goal is to assess/detect if two distant regions are linked by more read pairs than expected (in the controls).
##' @title Pair discordant reads
##' @param chr the chromosome name
##' @param start.pos the start position
##' @param end.pos the end position
##' @param files.df a data.frame with the path to the bam files for each samples.
##' @param samp the name of sample of interest.
##' @param controls the names of the samples to use as control.
##' @param bins.df a data.frame with the bins to use. If NULL (default), 1 Kbp bins will be used.
##' @param plot should a graph be displayed. Default is TRUE.
##' @param nb.cores the number of cores to use. Default is 1.
##' @return a data.frame with the mapping summary for each pair of region.
##' @author Jean Monlong
##' @export
pair.discordant.reads <- function(chr, start.pos, end.pos, files.df, samp, controls, bins.df=NULL, plot=TRUE, nb.cores=1){

  if(!is.null(bins.df)){
    bins.df = fragment.genome.hg19(1e3)
  }

  if(!all(c(samp, controls) %in% files.df$sample)){
    stop("'sample' and 'controls' samples must all be in 'files.df'.")
  }

  gr = GenomicRanges::GRanges(chr, IRanges::IRanges(start.pos, end.pos))
  bins.df$bin = paste(bins.df$chr, bins.df$start, sep="-")
  bins.gr = with(bins.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
  
  pairBins.f <- function(samp){
    param = Rsamtools::ScanBamParam(which=gr,
                         what = c("mrnm","mpos"),
                         flag = Rsamtools::scanBamFlag(isPaired=TRUE,isProperPair=FALSE,isDuplicate=FALSE,isNotPassingQualityControls=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE)
                         )
    bam.file = files.df$bam[which(files.df$sample==samp)]
    bam = Rsamtools::scanBam(bam.file,param=param)
    reads.gr = GenomicRanges::GRanges(bam[[1]]$mrnm, IRanges::IRanges(bam[[1]]$mpos, width=1))
    pairs.df = bins.df
    pairs.df$read = GenomicRanges::countOverlaps(bins.gr, reads.gr)
    pairs.df$sample = samp
    return(pairs.df)
  }

  pairs.df = as.data.frame(data.table::rbindlist(parallel::mclapply(c(samp, controls), pairBins.f, mc.cores=nb.cores)))

  ## Remove non-covered bins
  covered.bins = unique(pairs.df$bin[which(pairs.df$read>0)])
  pairs.df = pairs.df[which(pairs.df$bin %in% covered.bins), ]

  if(plot){
    bin = read = NULL ## Uglily appeases R checks
    print(ggplot2::ggplot(pairs.df[which(pairs.df$sample!=samp),], ggplot2::aes(x=bin, y=read, group=bin)) +
          ggplot2::geom_boxplot() + ggplot2::geom_point(data=pairs.df[which(pairs.df$sample==samp),], size=2, colour="blue") +
          ggplot2::theme_bw() + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1), legend.position="bottom") +
          ggplot2::xlab("genomic region") + ggplot2::ggtitle(paste("Discordant read-pairs linking to",chr,":",start.pos,"-",end.pos)))
  }
    
  return(pairs.df)
}
