##' Plot the raw read coverage in a genomic regions across samples.
##' @title  Plot the raw coverage in a region
##' @param gr a 'GRanges' object corresponding to the region of interest
##' @param files.df the data.frame with the path to the bam files for each sample
##' @param samples the samples to represent
##' @param proper should proper mapping be counted. Default is TRUE.
##' @param map.quality the minimum mapping quality. Default is 30.
##' @param bp.res the resolution, in bp. Default is 1 bp.
##' @return a ggplot object
##' @author Jean Monlong
##' @export
coverage.plot.raw <- function(gr, files.df, samples, proper=TRUE, map.quality=30, bp.res=1, flanks=NULL){

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
  cov.gr <- function(gr.i,bin.size=10){
    gr.bk = seq(min(start(gr.i)), max(end(gr.i)), bin.size)
    gr.cov = GRanges(seqnames(gr.i)[1], IRanges(gr.bk[-length(gr.bk)], width=bin.size))
    cov.t = table(cut((start(gr.i)+end(gr.i))/2, gr.bk, labels=FALSE))
    gr.cov$cov = 0
    gr.cov$cov[as.integer(names(cov.t))] = as.integer(cov.t)
    gr.cov
  }

  if(!is.null(flanks)){
    bkpts = c(start(gr), end(gr))
    start(gr) = start(gr)-flanks
    end(gr) = end(gr)+flanks
  }

  if(bp.res>1){
    gr.l = lapply(samples, function(samp.i){
                    cov.gr = cov.gr(gr.bam(files.df$bam[files.df$sample==samp.i], gr, proper=proper, map.quality=map.quality), bp.res)
                    cov.gr$sample = samp.i
                    cov.gr
                  })
    gr.ul = GRanges()
    for(ii in 1:length(gr.l)) gr.ul = c(gr.ul, gr.l[[ii]])
    final.p = ggbio::autoplot(gr.ul, geom="area",ggplot2::aes(y=cov),facets=sample~.) + ggplot2::theme_bw() + ggplot2::facet_grid(sample~.,scales="free")
  } else {
    gr.l = lapply(samples, function(samp.i){
                    reads.gr = gr.bam(files.df$bam[files.df$sample==samp.i], gr, proper=proper, map.quality=map.quality)
                    reads.gr$sample = samp.i
                    reads.gr
                  })
    ##gr.l = unlist(GenomicRanges::GRangesList(gr.l))
    gr.ul = GRanges()
    for(ii in 1:length(gr.l)) gr.ul = c(gr.ul, gr.l[[ii]])
    gr.ul$sample = factor(gr.ul$sample, levels=samples)
    final.p = ggbio::autoplot(gr.ul, facets=sample~.,stat="coverage") + ggplot2::theme_bw() + ggplot2::facet_grid(sample~.,scales="free")
  }

  if(!is.null(flanks)){
    final.p = final.p + ggplot2::geom_vline(xintercept=bkpts,linetype=2)
  }

  return(final.p)
}
