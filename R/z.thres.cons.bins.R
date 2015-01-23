##' Find a threshold on the Z-scores associated with how many consecutive bins are found.
##' @title Z-score thresholding using bin consecutiveness.
##' @param z.df a data.frame with at least 'chr', 'start', 'end' and 'z' columns.
##' @param plot should some graphs be displayed. Default if FALSE.
##' @param pvalues should the P-values be deduced from these thresholds. Default is FALSE. Note: the P-values and Q-values derive from custom 'recipes' and should be taken as informative only, for further filtering or priorization. 
##' @return a list with:
##' \item{dup.th}{the threshold for positive Z-scores (duplication signal)}
##' \item{del.th}{the threshold for negative Z-scores (deletion signal)}
##' \item{nb.ab.bins}{the number of abnormal bins}
##' \item{prop.ab.bins}{the proportion of abnormal bins}
##' \item{z.df}{a data.frame, similar to the input 'z.df' with new columns: 'abnormal' for bins with abnormal Z-scores and potentially 'pv'/'qv' for P-values/Q-values.}
##' @author Jean Monlong
##' @keywords internal
z.thres.cons.bins <- function(z.df, plot=FALSE, pvalues=FALSE){

  bin.w = round(median(z.df$end-z.df$start, na.rm=TRUE))
  cons.dist.f <- function(df){
    gr = with(df, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    gr.r = GenomicRanges::reduce(gr, min.gapwidth=2)
    data.frame(nbc = round(GenomicRanges::width(gr.r)/bin.w))
  }
  localMax <- function(x,y=NULL,min.max.prop=.1, loc.max=TRUE){
    if(is.null(y)){
      d = density(x,na.rm=TRUE)
    } else {
      d = data.frame(x=x, y=y)
      d = dplyr::arrange(d, x)
    }
    my = max(d$y)
    if(!loc.max){
      d$y = my-d$y
      my = max(d$y)      
    }
    im = 1+which(diff(sign(diff(d$y)))==-2)
    max.id = im[which(d$y[im] >= min.max.prop * my)]
    if(length(max.id)==0){
      return(list(lM=NA, h=NA))
    }
    max.id.o = max.id[order(d$y[max.id],decreasing=TRUE)]
    return(list(lM=d$x[max.id.o], h=d$y[max.id.o]/my))
  }
  find.th <- function(df, z.int=seq(1,20,.2)){
    nbcc.df = plyr::ldply(z.int,function(z.th){
      df =  dplyr::do(dplyr::group_by(subset(df, abs(z)>z.th), chr), cons.dist.f(.))
      df =  dplyr::summarize(dplyr::group_by(df, nbc), n=n())
      df$z.th=z.th
      df$p=df$n/sum(df$n)
      df
    })
    z.th = c(min(localMax(subset(nbcc.df, nbc==1)$z.th, subset(nbcc.df, nbc==1)$p)$lM),
      sapply(2:3, function(nbc.i)min(localMax(subset(nbcc.df, nbc==nbc.i)$z.th, subset(nbcc.df, nbc==nbc.i)$p, loc.max=FALSE)$lM)))
    max(z.th, na.rm=TRUE)
  }

  ## Split between duplication/deletion signal
  dup.df = subset(z.df, z>0)
  del.df = subset(z.df, z<0)

  ## Find threshold; second run scan with more resolution.
  dup.th = find.th(dup.df)
  ##dup.th = find.th(dup.df, seq(dup.th-.5, dup.th+.5, .01))
  ##dup.th = find.th(dup.df, seq(dup.th-.1, dup.th+.1, .005))
  del.th = find.th(del.df)
  ##del.th = find.th(del.df, seq(del.th-.5, del.th+.5, .01))
  ##del.th = find.th(del.df, seq(del.th-.1, del.th+.1, .005))

  ## P-value computation
  if(pvalues){
    sd.int = seq(0,dup.th,.01)
    pv.sd.dup = sapply(sd.int, function(sd.i) pnorm(dup.th, sd=sd.i, lower.tail = FALSE))
    sd.dup = sd.int[which.min(abs(.005-pv.sd.dup))]
    dup.df$pv = 2*pnorm(-dup.df$z,0,sd.dup)
    sd.int = seq(0,del.th,.01)
    pv.sd.del = sapply(sd.int, function(sd.i) pnorm(del.th, sd=sd.i, lower.tail = FALSE))
    sd.del = sd.int[which.min(abs(.005-pv.sd.del))]
    del.df$pv = 2*pnorm(del.df$z,0,sd.del)
  } 

  ## Annotate data.frames
  dup.df$abnormal = dup.df$z > dup.th
  del.df$abnormal = del.df$z < -del.th
  z.df = rbind(dup.df, del.df)
  
  ## Some statistics on the calls
  nb.ab.bins = sum(z.df$abnormal)
  prop.ab.bins = mean(z.df$abnormal)

  ## Multiple test correction
  if(pvalues){
    if(any(z.df$pv==0))
      z.df$pv[z.df$pv==0] = .Machine$double.xmin
    z.df$qv = p.adjust(z.df$pv,method="fdr")
  }

  ## Z-score distribution with thresholds
  if(plot){
    print(ggplot2::ggplot(z.df, ggplot2::aes(x=z)) + ggplot2::geom_histogram() + ggplot2::xlim(-20,20) + ggplot2::geom_vline(xintercept=c(-del.th, dup.th), linetype=2) + ggplot2::theme_bw())
    if(pvalues){
      print(ggplot2::ggplot(z.df, ggplot2::aes(x=pv)) + ggplot2::geom_histogram() + ggplot2::xlim(0,1) + ggplot2::theme_bw())
      print(ggplot2::ggplot(z.df, ggplot2::aes(x=qv)) + ggplot2::geom_histogram() + ggplot2::xlim(0,1) + ggplot2::theme_bw())
    }
  }
  
  return(list(dup.th=dup.th, del.th=del.th,nb.ab.bins=nb.ab.bins, prop.ab.bins=prop.ab.bins,
              z.df=subset(z.df,abnormal)))
}
