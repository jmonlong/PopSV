##' All the ranges defined by columns 'chr', 'start' and 'end' in the input data.frame are overlapped and used to define sub-ranges. Then the number of input ranges overlapping each sub-ranges is computed and return in a data.frame.
##' @title Frequency computation for ranges
##' @param range.df a data.frame with columns 'chr', 'start' and 'end'.
##' @param plot should a graph with the frequency distribution be displayed.
##' @param annotate.only If TRUE, the input ranges are annotated with the number of other overlapping ranges. If FALSE (Default), the input ranges are fragmented and the frequency computed for each sub-range.
##' @return a data.frame with the frequency of each sub-range
##' @author Jean Monlong
##' @export
##' @import magrittr
freq.range <- function(range.df, plot=FALSE, annotate.only=FALSE){
  if(!all(c("chr","start","end") %in% colnames(range.df))){
    stop("Missing column in 'range.df'. 'chr', 'start' and 'end' are required.")
  }
  nb.samp = length(unique(range.df$sample))
  freq.chr.gr <- function(cnv.o){
    gr =  with(cnv.o, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    if(annotate.only) {
      gr.d = gr
    } else {
      gr.d = GenomicRanges::disjoin(gr)
      cnv.o = GenomicRanges::as.data.frame(gr.d)[,1:3]
      colnames(cnv.o)[1] = "chr"
    }
    ol = GenomicRanges::findOverlaps(gr.d, gr)
    thits = base::table(IRanges::queryHits(ol))
    cnv.o$nb = 0
    cnv.o$nb[as.numeric(names(thits))] = as.numeric(thits)
    cnv.o$prop = cnv.o$nb / nb.samp
    cnv.o
  }
  fr.df = freq.chr.gr(range.df)
  if(plot){
    chr = nb = prop = gen.kb = NULL ## Uglily silence R checks
    f.df = fr.df %>% dplyr::group_by(chr, start, end) %>% dplyr::summarize(nb=sum(nb), prop=sum(prop), gen.kb=head((end-start)/1e3, 1)) %>% dplyr::arrange(chr)
    print(suppressWarnings(ggplot2::ggplot(f.df, ggplot2::aes(x=signif(prop,3), y=gen.kb, fill=chr)) + ggplot2::xlab("proportion of samples") + 
          ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() +
          ggplot2::ylab("abnormal genome (Kb)") +
          ggplot2::guides(fill=FALSE) + ggplot2::scale_x_continuous(breaks=seq(0,1,.2)) + 
          ggplot2::scale_fill_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),3))
          ))
  }
  return(fr.df)
}
