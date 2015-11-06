##' Plot the copy number estimates in a genomic regions across samples.
##' @title  Plot a region copy-number estimates
##' @param chr the chromosome
##' @param start the start position of the region
##' @param end the end position of the region
##' @param files.df a data.frame with the path to the different files associated with each sample. If 'sample' is absent from 'bc.f' the correct file from 'files.df' will be used. Default is NULL.
##' @param samples the set of samples to represent.
##' @param flanks the size of the flanking region. Default is 10000.
##' @param highlight.region Should the input region highlighted. Default is TRUE.
##' @return a ggplot object
##' @author Jean Monlong
##' @export
cn.plot <- function(chr, start, end, files.df, samples=NULL, flanks=1e4, highlight.region=TRUE){
  if(!is.null(samples)){
    files.df = subset(files.df, sample %in% samples)
  }
  exist = which(file.exists(files.df$fc) | file.exists(files.df$fc.gz))
  samples = files.df$sample[exist]
  if(length(samples)==0) {
    stop("No suitable samples found. Either not present in 'files.df' or 'fc' files not found.")
  }

  gr = GenomicRanges::GRanges(chr, IRanges::IRanges(start-flanks, end+flanks))

  cn.l = lapply(samples, function(samp){
                  cn = read.bedix(files.df$fc[which(files.df$sample==samp)], gr)
                  cn$sample = samp
                  colnames(cn)[4] = "value"
                  cn$pos = as.numeric(with(cn, (start+end)/2))
                  cn[,c("chr","sample","value","pos")]
                })
  cn.df = do.call(rbind, cn.l)

  max.cn = max(cn.df$value, na.rm=TRUE)
  pos = value = ggpSck = NULL ## Uglily appease R checks
  gp.o = ggplot2::ggplot(cn.df)
  if(highlight.region){
    gp.o = gp.o + ggplot2::geom_rect(xmin=start, xmax=end, ymin=0, ymax=max.cn, fill="yellow2", ggplot2::aes(alpha=ggpSck), data=data.frame(ggpSck=0)) + ggplot2::guides(alpha=FALSE)
  }
  gp.o = gp.o + ggplot2::theme_bw() + ggplot2::ylab("estimated copy number") + ggplot2::xlab("position") + ggplot2::geom_violin(ggplot2::aes(x=pos, y=2*value, group=pos), scale="width")

  gp.o
}
