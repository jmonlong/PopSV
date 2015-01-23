##' Chromosomes with clear aneuploidy are detected using a simple threshold-based approach. It is used to flag chromosome with complete aneuploidy that could handicap later analysis. It might not detect partial chromosomal aberration or aneuploidy in samples with noisy read coverage.
##' @title Flag chromosomal aneuploidy
##' @param samp the name of the sample to analyze.
##' @param files.df a data.frame with information about samples such as the path to the count files.
##' @param col.file the name of the column in 'files.df' with the path information to use.
##' @param nb.bins the number of bins to use when subsampling each chromosome. Default is 1000.
##' @param prop.aneu the proportion of normal bins under which a chromosome is considered aneuploid.
##' @param plot should the read coverage per chromosome be displayed. Default is FALSE.
##' @return a vector with the names of the flagged(aneuploid) chromosomes.
##' @author Jean Monlong
##' @export
aneuploidy.flag <- function(samp, files.df, col.file="bc.gz", nb.bins=1e3, prop.aneu=.1, plot=FALSE){
  
  localMax <- function(x,min.max.prop=.1,max=FALSE){
    d = density(x,na.rm=TRUE)
    im = 1+which(diff(sign(diff(d$y)))==-2)
    my = max(d$y)
    max.id = im[which(d$y[im] >= min.max.prop * my)]
    max.id.o = max.id[order(d$y[max.id],decreasing=TRUE)]
    return(list(lM=d$x[max.id], h=d$y[max.id]/my))
  }

  df = read.table(files.df[files.df$sample==samp,col.file], header=TRUE, as.is=TRUE)
  chrs = unique(df$chr)
  df = subset(df, bc>0)
  
  df.sub = dplyr::do(dplyr::group_by(df, chr), {.[sample.int(nrow(.),nb.bins),]})
  lm.o = localMax(df.sub$bc)
  lm.o = lm.o$lM[which.max(lm.o$h)]
  core.chrs = df.sub$chr[order(abs(lm.o-df.sub$bc))[1:(nrow(df.sub)*.3)]]
  cchrs.t = table(core.chrs)
  aneu.chrs = NULL
  if(length(cchrs.t)<length(chrs)){
    aneu.chrs = setdiff(chrs, names(cchrs.t))
  }
  if(any(cchrs.t < prop.aneu*nb.bins*.3)){
    aneu.chrs = c(aneu.chrs, names(cchrs.t)[which(cchrs.t < prop.aneu*nb.bins*.3)])
  }

  if(plot){
    df$aneu.flag = df$chr %in% aneu.chrs
    df$bc = winsor(df$bc, u=quantile(df$bc, probs=.99, na.rm=TRUE))
    p1 = ggplot2::ggplot(df, ggplot2::aes(x=bc, fill=aneu.flag)) + ggplot2::geom_density(alpha=.4) + ggplot2::theme_bw() + ggplot2::facet_grid(chr~., scales="free") + ggplot2::xlab("raw read coverage") + ggplot2::ylab("bin density") + ggplot2::theme(axis.text.y=ggplot2::element_blank()) + ggplot2::theme(legend.position="bottom")
    df$all = "all"
    p2 =  ggplot2::ggplot(df, ggplot2::aes(x=bc)) + ggplot2::geom_density(fill="grey70", alpha=.4) + ggplot2::theme_bw() + ggplot2::ylab("bin density") + ggplot2::xlab("raw read coverage") + ggplot2::theme(axis.text.y=ggplot2::element_blank())+ ggplot2::facet_grid(all~., scales="free") + ggplot2::ggtitle(samp)
    grid::grid.newpage()
    print(p1, vp=grid::viewport(1, 0.8, x=0.5, y=0.4))
    print(p2, vp=grid::viewport(1, 0.2, x=0.5, y=0.9))
  }

  return(aneu.chrs)
}
