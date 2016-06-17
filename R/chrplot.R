##' Display a chromosome view of the CNVs across samples
##' @title Chromosom plot of CNVs
##' @param cnv.df a data.frame with the CNVs (columns 'chr', 'start', 'end') for each sample (column 'sample')
##' @param type the type of graph to display. 'stacked' (default) shows the frequency line-graph across the chromosome. 'sample' shows a heatmap with each row representing one sample.
##' @param bin.size the bin size. Default is 1e6. 
##' @param chr.df a data.frame with the chromosome boundaries.
##' @return a ggplot graph object
##' @author Jean Monlong
##' @export
chrplot <- function(cnv.df, type=c("sample", "stacked"), bin.size=5e5, chr.df=NULL){
  ## Uglily appeases R checks
  cnv = seg = . = chr = NULL
  
  bins.df = fragment.genome.hg19(bin.size, XY.chr=TRUE, quiet=TRUE)

  if(is.null(chr.df)){
    chr.df = lapply(unique(cnv.df$chr), function(chr.a){
      data.frame(chr=chr.a, start=min(bins.df$start[which(bins.df$chr==chr.a)]),
                 end=max(bins.df$end[which(bins.df$chr==chr.a)]))
    })
    chr.df = do.call(rbind, chr.df)
  }

  mergeDup <- function(df){
    cnv.rle = rle(df$cnv)
    df$seg = rep(1:length(cnv.rle$lengths), cnv.rle$lengths)
    df = dplyr::group_by(df, seg, cnv)
    df = dplyr::summarize(df, start=min(start), end=max(end))
    df = dplyr::ungroup(df)
    df$seg = NULL
    df
  }
  
  if(type=="sample"){
    binF <- function(df, bins.df){
      bins.df = bins.df[which(bins.df$chr %in% unique(df$chr)),]
      bins.df$cnv = IRanges::overlapsAny(GenomicRanges::makeGRangesFromDataFrame(bins.df), GenomicRanges::makeGRangesFromDataFrame(df))
      bins.df = dplyr::group_by(bins.df, chr)
      bins.df = dplyr::arrange(bins.df, start)
      bins.df = dplyr::do(bins.df, {mergeDup(.)})
      bins.df
    }
  }

  cnv.samp = dplyr::group_by(cnv.df, sample)
  cnv.b = dplyr::do(cnv.samp, {binF(., bins.df)})
  cnv.b$sample = factor(cnv.b$sample, levels=unique(cnv.b$sample))
  cnv.b = cnv.b[which(cnv.b$cnv),]

  sampChr = unique(cnv.b[, c("sample","chr")])
  chr.df = merge(chr.df, sampChr)

  ggplot2::ggplot(cnv.b, ggplot2::aes(xmin=start/1e6, xmax=end/1e6, ymin=as.numeric(sample)-.5, ymax=as.numeric(sample)+.5)) + ggplot2::geom_rect(data=chr.df, alpha=0, colour="red") + ggplot2::geom_rect() + ggplot2::scale_y_continuous(breaks=1:nlevels(cnv.b$sample), labels=levels(cnv.b$sample)) + ggplot2::theme_bw() + ggplot2::xlab("position (Mb)") + ggplot2::facet_grid(chr~., scales="free", labeller=ggplot2::label_both)
  
}
