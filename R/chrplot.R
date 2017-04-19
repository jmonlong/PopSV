##' Display a chromosome view of the CNVs across samples
##' @title Chromosom plot of CNVs
##' @param cnv.df a data.frame with the CNVs (columns 'chr', 'start', 'end') for each sample (column 'sample')
##' @param type the type of graph to display. 'stacked' (default) shows the frequency line-graph across the chromosome. 'sample' shows a heatmap with each row representing one sample.
##' @param bin.size the bin size. Default is 1e6.
##' @param chr.df a data.frame with the chromosome boundaries.
##' @param showSampleNames should the sample names be displayed. Default is FALSE (it quickly becomes unreadable).
##' @param chr.prefix Is there a 'chr' prefix in the chromosome names. Default is FALSE
##' @param show.all.chr Should all the chromosomes be displayed (even if no CNVs) ? Default is FALSE, i.e. only chr with CNVs are displayed.
##' @return a ggplot graph object
##' @author Jean Monlong
##' @export
chrplot <- function(cnv.df, type=c("sample", "stacked"), bin.size=5e5, chr.df=NULL, showSampleNames=FALSE, chr.prefix=FALSE, show.all.chr=FALSE){
  ## Uglily appeases R checks
  cnv = seg = . = chr = NULL

  bins.df = fragment.genome.hg19(bin.size, XY.chr=TRUE, quiet=TRUE, chr.prefix=chr.prefix)
  if(chr.prefix){
      bins.df$chr = factor(bins.df$chr, levels=paste0("chr",c(1:22,"X","Y")))
  } else {
      bins.df$chr = factor(bins.df$chr, levels=c(1:22,"X","Y"))
  }

  if(is.null(chr.df)){
    if(show.all.chr){
      chrs = unique(bins.df$chr)
    } else {
      chrs = unique(cnv.df$chr)
    }
    chr.df = lapply(chrs, function(chr.a){
      data.frame(chr=chr.a, start=min(bins.df$start[which(bins.df$chr==chr.a)]),
                 end=max(bins.df$end[which(bins.df$chr==chr.a)]))
    })
    chr.df = do.call(rbind, chr.df)
    if(chr.prefix){
        chr.df$chr = factor(chr.df$chr, levels=paste0("chr",c(1:22,"X","Y")))
    } else {
        chr.df$chr = factor(chr.df$chr, levels=c(1:22,"X","Y"))
    }
  }

  olProp <- function(qgr, sgr){
    sgr = reduce(sgr)
    ol = GenomicRanges::findOverlaps(qgr, sgr)
    cov.t = tapply(width(GenomicRanges::pintersect(sgr[S4Vectors::subjectHits(ol)],qgr[S4Vectors::queryHits(ol)])), S4Vectors::queryHits(ol), sum)
    out = rep(0,length(qgr))
    out[as.numeric(names(cov.t))] = as.numeric(cov.t)
    out / GenomicRanges::width(qgr)
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

  makeGRanges <- function(df){
    if(all(c("chr","start","end") %in% colnames(df))){
      df = df[,c("chr","start","end")]
    }
    GenomicRanges::makeGRangesFromDataFrame(df)
  }

  if(type[1]=="sample"){
    binF <- function(df, bins.df){
      bins.df = bins.df[which(bins.df$chr %in% unique(df$chr)),]
      ## bins.df$cnv = IRanges::overlapsAny(GenomicRanges::makeGRangesFromDataFrame(bins.df), GenomicRanges::makeGRangesFromDataFrame(df))
      bins.df$cnv = olProp(GenomicRanges::makeGRangesFromDataFrame(bins.df), GenomicRanges::makeGRangesFromDataFrame(df))
      bins.df = dplyr::group_by(bins.df, chr)
      bins.df = dplyr::arrange(bins.df, start)
      bins.df = dplyr::do(bins.df, {mergeDup(.)})
      bins.df
    }
  }

  cnv.samp = dplyr::group_by(cnv.df, sample)
  cnv.b = dplyr::do(cnv.samp, {binF(., bins.df)})
  cnv.b$sample = factor(cnv.b$sample, levels=unique(cnv.b$sample))
  cnv.b = cnv.b[which(cnv.b$cnv>0),]

  ## sampChr = unique(cnv.b[, c("sample","chr")])
  sampChr = data.frame(sample=rep(unique(cnv.b$sample),each=length(chrs)),
                       chr=rep(chrs, length(unique(cnv.b$sample))),
                       stringsAsFactors=FALSE)
  chr.df = merge(chr.df, sampChr)

  ggp = ggplot2::ggplot(cnv.b, ggplot2::aes(xmin=start/1e6, xmax=end/1e6, ymin=as.numeric(sample)-.5, ymax=as.numeric(sample)+.5)) + ggplot2::geom_rect(data=chr.df, alpha=.1, fill="black") + ggplot2::geom_rect(ggplot2::aes(fill=cnv)) + ggplot2::scale_y_continuous(breaks=1:nlevels(cnv.b$sample), labels=levels(cnv.b$sample)) + ggplot2::theme_bw() + ggplot2::xlab("position (Mb)") + ggplot2::facet_grid(chr~., scales="free") + ggplot2::scale_fill_gradient(name=paste(bin.size,"bp bin\noverlap"), high="red",low="blue", limits=0:1)

  if(!showSampleNames){
    ggp = ggp + ggplot2::theme(axis.text.y=ggplot2::element_blank())
  }

  if(length(unique(cnv.b$sample))==1){
    ggp = ggp + ggplot2::ggtitle(cnv.b$sample[1])
  }

  ggp
}
