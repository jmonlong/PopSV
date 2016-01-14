##' Summary of regions with abnormal coverage as detected by 'call.abnomal.cov'. It's a static version of 'sv.summary.interactive'. It produces graphs to to get an idea of the quality of the calls. For normal samples we would expect: the amount of called genome in each sample to be somewhat similar; the copy number estimates to locate near integer values; No systematic calls (i.e. called in all samples). The output graphs represent the amount of affected genome across samples, the copy number estimates distribution, the frequency distribution.
##'
##' A list is returned. Its first element is a list of ggplot object with the different graphs. The second is a data.frame with summary information used to create the graphs.
##' @title Abnormal regions summary
##' @param res.df the data.frame with the results.
##' @param out.pdf the name for the output PDF file to create. If NULL (default), the graphs are either displayed directly or returned as a list.
##' @param print Should the graphs be displayed ? If 'out.pdf' though, this parameter is not used as no graph will be displayed.
##' @return A list with :
##' \item{graphs.l}{a list of ggplot graphs}
##' \item{info.df}{a data.frame with different summary statistics}
##' @author Jean Monlong
##' @export
sv.summary <- function(res.df, out.pdf=NULL, print=TRUE){
  sample = gen.kb = col = cn = chr = nb = fc = type = start = end = prop = . = freq.n = nb.bin.cons = NULL ## Uglily appease R checks

  nb.samp = length(unique(res.df$sample))
  freq.chr.gr <- function(cnv.o){
    gr =  with(cnv.o, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    gr.d = GenomicRanges::disjoin(gr)
    ol = GenomicRanges::findOverlaps(gr.d, gr)
    thits = base::table(IRanges::queryHits(ol))
    freq.df = data.frame(win.id = as.numeric(names(thits)), nb=as.numeric(thits))
    win.df = GenomicRanges::as.data.frame(gr.d[as.numeric(freq.df$win.id)])[,1:3]
    freq.df$chr = as.character(win.df$seqnames)
    freq.df$start = as.numeric(win.df$start)
    freq.df$end = as.numeric(win.df$end)
    freq.df$prop = freq.df$nb / nb.samp
    freq.df
  }

  res.df$gen.kb = with(res.df, (end-start)/1e3)

  samp.o = aggregate(gen.kb~sample, data=res.df, sum)
  res.df$sample = factor(res.df$sample, levels=samp.o$sample[order(samp.o$gen.kb)])

  info.df = data.frame(variable="kb-meanPerSamp",value=mean(samp.o$gen.kb), stringsAsFactors=FALSE)

  freq.df = rbind(data.frame(type="deletion",freq.chr.gr(res.df[which(res.df$fc<1),])),
          data.frame(type="duplication",freq.chr.gr(res.df[which(res.df$fc>1),])))

  output = list()

  res.df$col = ifelse(res.df$fc>1, "duplication","deletion")
  res.df.s = dplyr::summarize(dplyr::group_by(res.df, sample, col), gen.kb=sum(gen.kb))
  output$nb.calls.perType = ggplot2::ggplot(res.df.s, ggplot2::aes(x=sample, y=gen.kb, fill=col)) +
    ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() + ggplot2::scale_fill_brewer(name="type",palette="Set1") +
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     legend.position=c(0,1),legend.justification=c(0,1)) + ggplot2::ylab("abnormal genome (Kb)")

  m.df = dplyr::summarize(dplyr::group_by(res.df.s,col), gen.kb=mean(gen.kb))
  info.df = rbind(info.df, data.frame(variable=paste("kb-meanPerSamp",m.df$col,sep="-"), value=m.df$gen.kb, stringsAsFactors=FALSE))

  res.df$col = cut(res.df$nb.bin.cons, breaks=c(0,1,2,3,5,10,Inf))
  res.df.s = dplyr::summarize(dplyr::group_by(res.df, sample, col), gen.kb=sum(gen.kb))
  output$nb.calls.perSize = ggplot2::ggplot(res.df.s, ggplot2::aes(x=sample, y=gen.kb, fill=col)) +
    ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() + ggplot2::scale_fill_brewer(name="size (bin)",palette="Set1") +
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     legend.position=c(0,1),legend.justification=c(0,1)) + ggplot2::ylab("abnormal genome (Kb)")

  pdf = res.df[which(res.df$nb.bin.cons>=3),]
  pdf$cn = pdf$fc*2
  if(any(pdf$cn>5)) pdf$cn[pdf$cn>5] = 5
  pdf$col = ifelse(pdf$fc>1, "duplication","deletion")
  output$cn.perType =  ggplot2::ggplot(pdf, ggplot2::aes(x=cn, fill=col)) +
    ggplot2::geom_histogram() + ggplot2::theme_bw() + ggplot2::scale_x_continuous(breaks=seq(0,5,1), labels=c(seq(0,4,1),">5"), limits=c(0,5.2)) + ggplot2::scale_fill_brewer(name="type",palette="Set1") +
      ggplot2::theme(legend.position=c(.5,1),legend.justification=c(.5,1)) + ggplot2::xlab("Copy Number estimates") + ggplot2::ylab("number of calls")

  pdf$col = cut(pdf$nb.bin.cons, breaks=c(0,1,2,3,5,10,Inf))
  output$cn.perSize =  ggplot2::ggplot(pdf, ggplot2::aes(x=cn, fill=col)) +
    ggplot2::geom_histogram() + ggplot2::theme_bw() + ggplot2::scale_x_continuous(breaks=seq(0,5,1), labels=c(seq(0,4,1),">5"), limits=c(0,5.2)) + ggplot2::scale_fill_brewer(name="size (bin)",palette="Set1") +
      ggplot2::theme(legend.position=c(.5,1),legend.justification=c(.5,1)) + ggplot2::xlab("Copy Number estimates") + ggplot2::ylab("number of calls")


  f.df = dplyr::summarize(dplyr::group_by(freq.df, chr, start, end), nb=sum(nb), prop=sum(prop), gen.kb=head((end-start)/1e3, 1))
  f.df = dplyr::summarize(dplyr::group_by(f.df, chr, nb, prop), gen.kb=sum(gen.kb))
  output$freq = ggplot2::ggplot(dplyr::arrange(f.df, chr), ggplot2::aes(x=nb,  y=gen.kb, fill=chr)) + ggplot2::xlab("number of samples") +
    ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() + ggplot2::ylab("abnormal genome (Kb)") +
      ggplot2::guides(fill=FALSE) + ggplot2::scale_fill_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),3))

  output$freq.10samps = ggplot2::ggplot(dplyr::arrange(f.df[which(f.df$nb>=10),], chr), ggplot2::aes(x=nb,  y=gen.kb, fill=chr)) + ggplot2::xlab("number of samples") +
    ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() + ggplot2::ylab("abnormal genome (Kb)") +
      ggplot2::guides(fill=FALSE) + ggplot2::scale_fill_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),3))

  info.df = rbind(info.df, data.frame(variable="freq<10%-propBp", value=sum(f.df$gen.kb[which(f.df$prop<.1)])/sum(f.df$gen.kb), stringsAsFactors=FALSE))

  if(!is.null(out.pdf)){
    pdf(out.pdf, 9,6)
    lapply(output, print)
    dev.off()
  } else if(print) {
    par(ask=TRUE)
    lapply(output, print)
  }

  return(list(graphs.l = output, info.df=info.df))
}
