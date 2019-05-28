##' The frequency is defined as the number of input ranges overlapping a particular bp/region. In practice, all the ranges defined by columns 'chr', 'start' and 'end' in the input data.frame are overlapped and used to define sub-ranges. Then the number of input ranges overlapping each sub-ranges is computed and return in a data.frame.
##'
##' If 'annotate.only=TRUE' however, the frequency of an input range is computed as the number of ranges overlapping it. It is less-suited to describe the frequency distribution, but more convenient to filter out frequent variants.
##' @title Frequency computation for ranges
##' @param range.df a data.frame with columns 'chr', 'start' and 'end'.
##' @param plot should a graph with the frequency distribution be displayed. Default is FALSE.
##' @param annotate.only If TRUE, the input ranges are annotated with the number of other overlapping ranges. If FALSE (Default), the input ranges are fragmented and the frequency computed for each sub-range.
##' @param nb.samp the total number of sample. Used if not NA. Default is NA.
##' @param min.rol minimum reciprocal overlap when 'annotate.only=TRUE'. If 0 (default) any overlap counts.
##' @param chunk.size the number or regions in each chunk when 'annotate.only=TRUE'. Default is 1e5. Decrease if memory problems.
##' @return a data.frame with the frequency of each sub-range
##' @author Jean Monlong
##' @export
##' @import magrittr
freq.range <- function(range.df, plot=FALSE, annotate.only=FALSE, nb.samp=NA, min.rol=0, chunk.size=1e5){
  . = V1 = chr = nb = prop = gen.kb = queryHits = subjectHits = qw = sw = qsw = n = chunk = NULL ## Uglily silence R checks
  if(!all(c("chr","start","end") %in% colnames(range.df))){
    stop("Missing column in 'range.df'. 'chr', 'start' and 'end' are required.")
  }
  if(all(colnames(range.df)!="sample")){
    range.df$sample = 1:nrow(range.df)
  }
  if(is.na(nb.samp)){
    nb.samp = length(unique(range.df$sample))
  }
  gr.all = with(range.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), sample = sample))
  if(annotate.only){
    freq.chunk <- function(df) {
      gr = GenomicRanges::makeGRangesFromDataFrame(df)
      ol = GenomicRanges::findOverlaps(gr, gr.all)
      ol = GenomicRanges::as.data.frame(ol)
      ol = dplyr::mutate(ol, qw=GenomicRanges::width(gr)[queryHits],
          sw=GenomicRanges::width(gr.all)[subjectHits],
          qsw=GenomicRanges::width(GenomicRanges::pintersect(gr[queryHits], gr.all[subjectHits])),
          sample=gr.all$sample[subjectHits])
      ol = dplyr::filter(ol, qsw/qw>= min.rol, qsw/sw>= min.rol)
      ol = ol %>% dplyr::group_by(queryHits) %>% dplyr::summarize(nb=length(unique(sample)))
      df$nb = 0
      df$nb[ol$queryHits] = ol$nb
      df$prop = df$nb/nb.samp
      df
    }
    fr.df = range.df %>% dplyr::ungroup(.) %>%
      dplyr::mutate(chunk=cut(1:dplyr::n(), 1+ceiling(dplyr::n()/chunk.size))) %>% dplyr::group_by(chunk) %>%
      dplyr::do(freq.chunk(.)) %>% dplyr::ungroup(.) %>% dplyr::select(-chunk) %>% as.data.frame
  } else {
    gr.d = GenomicRanges::disjoin(gr.all)
    fr.df = GenomicRanges::as.data.frame(gr.d)[, 1:3]
    colnames(fr.df)[1] = "chr"
    ol = data.table::data.table(as.data.frame(GenomicRanges::findOverlaps(gr.d,gr.all)))
    ol = ol[, .(length(unique(gr.all$sample[subjectHits]))), by = .(queryHits)]
    fr.df$nb = 0
    fr.df$nb[ol[, queryHits]] = ol[, V1]
    fr.df$prop = fr.df$nb/nb.samp
  }
  if(plot){
    f.df = fr.df %>% dplyr::group_by(chr, start, end) %>% dplyr::summarize(nb=sum(nb), prop=sum(prop), gen.kb=utils::head((end-start)/1e3, 1)) %>% dplyr::arrange(chr)
    print(suppressWarnings(ggplot2::ggplot(f.df, ggplot2::aes(x=signif(prop,3), y=gen.kb, fill=chr)) + ggplot2::xlab("proportion of samples") +
                             ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() +
                               ggplot2::ylab("abnormal genome (Kb)") +
                                 ggplot2::guides(fill=FALSE) + ggplot2::scale_x_continuous(breaks=seq(0,1,.2)) +
                                   ggplot2::scale_fill_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),3))
                           ))
  }
  return(fr.df)
}
