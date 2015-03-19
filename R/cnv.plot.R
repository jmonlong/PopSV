##' Plot the coverage in a genomic regions across samples. The files ('bc.f', 'norm.stats.f') must be ordered and indexed.
##' @title  Plot a region coverage
##' @param chr the chromosome
##' @param start the start position of the region
##' @param end the end position of the region
##' @param bc.f the path to the bin count file.
##' @param norm.stats.f the path to the normalization statistics file.
##' @param sv.df the data.frame with the calls.
##' @param ref.samples the names of the reference samples.
##' @param boxplot should the reference be represented as boxplots. If FALSE, violin plots will be used.
##' @param samples the set of samples to represent.
##' @return a ggplot object
##' @author Jean Monlong
##' @export
cnv.plot <- function(chr, start, end, bc.f, norm.stats.f=NULL, sv.df=NULL, ref.samples=NULL, boxplot=FALSE , samples=NULL){

  gr = GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
  
  ## Read coverage file
  message("Bin count...")
  bc = read.bedix(bc.f, gr)
  bc.all = reshape::melt.data.frame(bc, id.vars=c("chr","start","end"),variable_name="sample")
  bc.all$pos = as.numeric(with(bc.all, (start+end)/2))
  
  if(is.null(ref.samples)){
    ref.samples = unique(bc.all$sample)
  }

  bc.sv = NULL
  if(!is.null(sv.df) & !is.null(nrow(sv.df))){
    message("PopSV calls...")
    bin.gr = unique(with(bc.all, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))))
    sv.gr =  with(sv.df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    sv.df = sv.df[overlapsAny(sv.gr, gr),]
    bc.sv = bc.all[which(bc.all$sample %in% unique(sv.df$sample)), ]
  }
  if(!is.null(samples)){
    bc.sv = bc.all[which(bc.all$sample %in% samples),]
  }
    
  if(!is.null(norm.stats.f)){
    bc.ref = read.bedix(norm.stats.f, gr)
    bc.ref$pos = as.numeric(with(bc.ref, (start+end)/2))
  } else {
    bc.ref = bc.all[which(bc.all$sample %in% ref.samples),]
  }

  ## Plot reference samples
  if(!is.null(norm.stats.f)) {
    gp.o = ggplot2::ggplot(bc.ref, ggplot2::aes(x=pos)) + ggplot2::theme_bw() + ggplot2::ylab("normalized coverage") + ggplot2::geom_errorbar(ggplot2::aes(ymin=m-3*sd,ymax=m+3*sd)) + ggplot2::geom_point(ggplot2::aes(y=m))
  } else if(boxplot){
    gp.o = ggplot2::ggplot(bc.ref, ggplot2::aes(x=pos, y=value)) + ggplot2::theme_bw() + ggplot2::ylab("normalized coverage") + ggplot2::geom_boxplot(ggplot2::aes(group=pos), fill="lightblue",alpha=.7)
  } else {
    gp.o = ggplot2::ggplot(bc.ref, ggplot2::aes(x=pos, y=value)) + ggplot2::theme_bw() + ggplot2::ylab("normalized coverage") + ggplot2::geom_violin(ggplot2::aes(group=pos), fill="lightblue",alpha=.7)
  }
  
  
  ## Add the SVs
  if(!is.null(bc.sv)){
    gp.o = gp.o + ggplot2::geom_line(ggplot2::aes(y=value, group=sample),alpha=.7,data=bc.sv)
  }

  gp.o  
}
