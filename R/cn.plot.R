##' Plot the copy number estimates in a genomic regions across samples.
##' @title  Plot copy-number estimates in a region.
##' @param files.df a data.frame with the path to the different files associated with each sample.
##' @param gr a GRange object with the region to represent. Default is NULL.
##' @param chr the chromosome. Used if 'gr' is NULL.
##' @param start the start position of the region. Used if 'gr' is NULL.
##' @param end the end position of the region. Used if 'gr' is NULL.
##' @param samples the set of samples to represent.
##' @param flanks the size of the flanking region. Default is 10000.
##' @param highlight.region Should the input region highlighted. Default is TRUE.
##' @param draw.lines Should lines be drawn to represent each sample. Default is TRUE.
##' @param nb.cores number of cores to use. Default is 1. Increase number if many samples are to be used.
##' @return a ggplot object
##' @author Jean Monlong
##' @export
cn.plot <- function(files.df, gr=NULL, chr=NULL, start=NULL, end=NULL, samples=NULL, flanks=1e4, highlight.region=TRUE, draw.lines=TRUE, nb.cores=1){
  if(!is.null(samples)){
    files.df = files.df[which(files.df$sample %in% samples),]
  }
  exist = which(file.exists(files.df$fc) | file.exists(files.df$fc.gz))
  samples = files.df$sample[exist]
  if(length(samples)==0) {
    stop("No suitable samples found. Either not present in 'files.df' or 'fc' files not found.")
  }

  if(is.null(gr)){
    if(is.null(chr) | is.null(start) | is.null(end)){
      stop("Either 'gr' or 'chr'/'start'/'end' arguments must be given.")
    }
    gr = GenomicRanges::GRanges(chr, IRanges::IRanges(start-flanks, end+flanks))
  } else {
    start = GenomicRanges::start(gr)[1]
    end = GenomicRanges::end(gr)[1]
    GenomicRanges::start(gr) = GenomicRanges::start(gr) - flanks
    GenomicRanges::end(gr) = GenomicRanges::end(gr) + flanks
  }

  cn.l = parallel::mclapply(samples, function(samp){
    cn = read.bedix(files.df$fc[which(files.df$sample==samp)], gr)
    cn$sample = samp
    colnames(cn)[4] = "value"
    cn$pos = as.numeric(with(cn, (start+end)/2))
    cn[,c("chr","sample","value","pos")]
  }, mc.cores=nb.cores)
  cn.df = do.call(rbind, cn.l)

  max.cn = ceiling(max(2*cn.df$value, na.rm=TRUE))
  pos = value = ggpSck = NULL ## Uglily appease R checks
  gp.o = ggplot2::ggplot(cn.df)
  if(highlight.region){
    gp.o = gp.o + ggplot2::geom_rect(xmin=start, xmax=end, ymin=0, ymax=max.cn, fill="yellow2", ggplot2::aes(alpha=ggpSck), data=data.frame(ggpSck=0)) + ggplot2::guides(alpha=FALSE)
  }
  if(draw.lines){
    gp.o = gp.o + ggplot2::geom_line(ggplot2::aes(x=pos, y=2*value, group=sample), alpha=.4)
  }
  if(max.cn<20){
    gp.o = gp.o + ggplot2::geom_hline(yintercept=0:max.cn, linetype=2)
  }
  gp.o = gp.o + ggplot2::theme_bw() + ggplot2::ylab("estimated copy number") + ggplot2::xlab("position") + ggplot2::geom_violin(ggplot2::aes(x=pos, y=2*value, group=pos), alpha=.6, scale="width") + ggplot2::ylim(0,max.cn)

  gp.o
}
