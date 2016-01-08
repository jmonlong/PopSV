##' Plot the bin coverage in a genomic regions across samples. The files ('bc.f', 'norm.stats.f') must be ordered and indexed.
##' @title  Plot the bin coverage in a region
##' @param chr the chromosome
##' @param start the start position of the region
##' @param end the end position of the region
##' @param bc.f the path to the bin count file.
##' @param norm.stats.f the path to the normalization statistics file.
##' @param sv.df the data.frame with the calls.
##' @param ref.samples the names of the reference samples.
##' @param boxplot should the reference be represented as boxplots. If FALSE, violin plots will be used.
##' @param samples the set of samples to represent.
##' @param files.df a data.frame with the path to the different files associated with each sample. If 'sample' is absent from 'bc.f' the correct file from 'files.df' will be used. Default is NULL.
##' @param anno.df a data.frame with additional information (e.g. gene annotation) to be added to the graph.
##' @param anno.col the name of 'anno.df' column to use to differentiate elements in the graph. E.g. 'geneName' to color the genes in the graph.
##' @param flanks the size of the flanking region. Default is 10000.
##' @param absolute.position should the bins be placed in absolute position. Default is TRUE. If FALSE, distant bins might be displayed next to each other if no bins are available in between (useful for targeted sequencing).
##' @return a ggplot object
##' @author Jean Monlong
##' @export
coverage.plot <- function(chr, start, end, bc.f, norm.stats.f=NULL, sv.df=NULL, ref.samples=NULL, boxplot=FALSE , samples=NULL, files.df=NULL, anno.df=NULL, anno.col="geneName", flanks=1e4, absolute.position=TRUE){

  gr = GenomicRanges::GRanges(chr, IRanges::IRanges(start-flanks, end+flanks))

  ## Read coverage file
  message("Bin count...")
  bc = read.bedix(bc.f, gr)
  bc.all = reshape::melt.data.frame(bc, id.vars=c("chr","start","end"),variable_name="sample")
  bc.all$pos = as.numeric(with(bc.all, (start+end)/2))
  start = min(bc.all$start[which(bc.all$end>start)])
  end = max(bc.all$end[which(bc.all$start<end)])

  if(is.null(ref.samples)){
    ref.samples = unique(bc.all$sample)
  }

  if(!is.null(norm.stats.f)){
    bc.ref = read.bedix(norm.stats.f, gr)
    bc.ref$pos = as.numeric(with(bc.ref, (start+end)/2))
    max.bc = max(bc.ref$m+3*bc.ref$sd, na.rm=TRUE)
  } else {
    bc.ref = bc.all[which(bc.all$sample %in% ref.samples),]
    max.bc = max(bc.ref$value, na.rm=TRUE)
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
    if(!all(samples %in% unique(bc.all$sample))){
      if(is.null(files.df)){
        warning("'files.df' is NULL and some samples are not in ",bc.f,", hence they will not be displayed.")
      } else {
        bc.sv.l = lapply(setdiff(samples, unique(bc.all$sample)), function(samp){
          bc = read.bedix(files.df$bc.gc.norm.gz[which(files.df$sample==samp)], gr)
          bc$sample = samp
          colnames(bc)[4] = "value"
          bc$pos = as.numeric(with(bc, (start+end)/2))
          bc[,c("chr","start","end","sample","value","pos")]
        })
        bc.sv = rbind(bc.sv, do.call(rbind, bc.sv.l))
      }
    }
  }
  if(!is.null(bc.sv)){
    max.bc = max(max.bc, max(bc.sv$value, na.rm=TRUE))
  }

  if(!absolute.position){
    bc.ref$pos = factor(round(bc.ref$pos))
  }
  
  ## Plot reference samples
  pos = m = value = ggpSck = NULL ## Uglily appease R checks
  gp.o = ggplot2::ggplot(bc.ref) + ggplot2::theme_bw() + ggplot2::ylab("normalized coverage") + ggplot2::xlab("position")
  if(flanks>0){
    gp.o = gp.o + ggplot2::geom_rect(xmin=start, xmax=end, ymin=0, ymax=max.bc, fill="yellow2", ggplot2::aes(alpha=ggpSck), data=data.frame(ggpSck=0)) + ggplot2::guides(alpha=FALSE)
  } 
  if(!is.null(norm.stats.f)) {
    gp.o = gp.o + ggplot2::geom_errorbar(ggplot2::aes(x=pos, ymin=m-3*sd,ymax=m+3*sd)) + ggplot2::geom_point(ggplot2::aes(y=m))
  } else if(boxplot){
    gp.o = gp.o + ggplot2::geom_boxplot(ggplot2::aes(x=pos, y=value, group=pos), fill="lightgreen",alpha=.7)
  } else {
    gp.o = gp.o + ggplot2::geom_violin(ggplot2::aes(x=pos, y=value, group=pos), fill="grey90",alpha=.7, scale="width")
  }

  ## Add the SVs
  if(!is.null(bc.sv)){
    gp.o = gp.o + ggplot2::geom_line(ggplot2::aes(x=pos, y=value, colour=sample, group=sample),alpha=.7, size=2, data=bc.sv) + ggplot2::geom_point(ggplot2::aes(x=pos, y=value, colour=sample),size=4, data=bc.sv)
  }

  if(!is.null(anno.df)){
    anno.df = anno.df[which(anno.df$chr==chr & anno.df$end>=start & anno.df$start<=end),]
    anno.df$y = factor(anno.df[,anno.col])
    anno.df$y = (as.numeric(anno.df$y)-1) * max.bc / nlevels(anno.df$y) /4
    gp.o = gp.o + ggplot2::geom_segment(ggplot2::aes_string(x="start",xend="end",y="y",yend="y",colour=anno.col), size=6, alpha=.5, data=anno.df)
  }

  if(!absolute.position){
    gp.o = gp.o + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
  }
  
  gp.o
}
