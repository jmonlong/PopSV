##' Plot the raw read coverage in a genomic regions across samples.
##' @title  Plot the raw coverage in a region
##' @param chr chromosome 
##' @param start genomic coordinate
##' @param end genomic coordinate
##' @param files.df a data.frame with the path to the bam files
##' @param samples the samples to represent
##' @param proper should proper mapping be counted. Default is TRUE.
##' @param map.quality the minimum mapping quality. Default is 30.
##' @param bp.res the resolution, in bp. Default is 1 bp.
##' @param flanks the amount of bp to add in the flanks.
##' @return a ggplot object
##' @author Jean Monlong
##' @export
coverage.plot.raw <- function(chr, start, end, files.df, samples, proper=TRUE, map.quality=30, bp.res=1, flanks=NULL){

  gr.bam <- function(bam.file, gr, proper = TRUE, map.quality = 30) {
    bai.file = sub("bam$", "bai", bam.file, perl = TRUE)
    if (!file.exists(bai.file)) {
      bai.file = paste0(bam.file, ".bai")
      if (!file.exists(bai.file)) {
        stop("Index file is missing (neither '.bai' nor '.bam.bai').")
      }
    }
    param = Rsamtools::ScanBamParam(which = gr, what = c("rname", "pos", "qwidth", "mapq"), flag = Rsamtools::scanBamFlag(isProperPair = proper, isDuplicate = FALSE, isNotPassingQualityControls = FALSE, isUnmappedQuery = FALSE))
    bam = Rsamtools::scanBam(bam.file, index = bai.file, param = param)
    bam = bam[which(unlist(lapply(bam, function(ee)length(ee[[1]])>0)))]
    bam.df = do.call(rbind, lapply(names(bam),function(bin) data.frame(bin=bin, as.data.frame(bam[[bin]]))))
    bam.df = bam.df[which(bam.df$mapq >= map.quality),]
    if(is.null(bam.df)) return(GenomicRanges::GRanges())
    return(with(bam.df, GenomicRanges::GRanges(rname, IRanges::IRanges(pos, width = qwidth), bin = bin)))
  }
  cov.reads <- function(reads.gr, win.gr){
    GenomicRanges::countOverlaps(win.gr, reads.gr)
  }
  
  if(!is.null(flanks)){
    bkpts = c(start, end)
    start = start-flanks
    end = end+flanks
  }

  gr.f = GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
  reads.l = lapply(samples, function(samp.i){
    gr.bam(files.df$bam[files.df$sample==samp.i], gr.f, proper=proper, map.quality=map.quality)
  })

  gr.bk = seq(GenomicRanges::start(gr.f), GenomicRanges::end(gr.f), bp.res)
  gr.frag = GenomicRanges::GRanges(chr, IRanges::IRanges(gr.bk[-length(gr.bk)], width = bp.res))
  cov.l = lapply(1:length(reads.l), function(ii) {
    cov.reads(reads.l[[ii]], gr.frag)
  })
  cov = matrix(0, length(gr.frag), length(samples))
  for (ii in 1:length(samples)) cov[, ii] = cov.l[[ii]]
  if(!is.null(flanks)){
    norm.f = apply(cov[c(1:floor(flanks/bp.res),(nrow(cov)-floor(flanks/bp.res)+1):nrow(cov)),],2,mean, na.rm=TRUE)
  } else {
    norm.f = rep(1,length(samples))
  }
  cov = cov %*% diag(mean(norm.f)/norm.f)
  colnames(cov) = samples
  rownames(cov) = GenomicRanges::start(gr.frag)
  cov.df = reshape::melt(cov)
  colnames(cov.df) = c("position","sample","cov")
  cov.df$sample = factor(cov.df$sample, levels=samples)
  cov.df$abnormal = cov.df$sample==samples[1]
  
  final.p = ggplot2::ggplot(subset(cov.df,!abnormal), ggplot2::aes(x=position,y=cov, group=sample)) + ggplot2::geom_line(alpha=.8) + ggplot2::geom_line(data=subset(cov.df,abnormal), size=3) + ggplot2::theme_bw() + ggplot2::scale_size_manual(values=c(1,2))
  if(!is.null(flanks)){
    final.p = final.p + ggplot2::geom_vline(xintercept=bkpts,linetype=2)
  }
  final.p
  
  return(final.p)
}
