##' Find breakpoints as the location with the most dramatic local coverage differences between a samples and a set of reference samples.
##' @title Breakpoint finder
##' @param bkpt.gr a GRanges object with one or several abnormal regions in a sample.
##' @param files.df a data.frame with information about the BAM file path for each sample (columns 'sample' and 'bam').
##' @param test.sample the name of the sample of interest.
##' @param ref.samples the names of the samples to use as reference for the breakpoint detection.
##' @param flank.size the size of the flanking region investigated on each side of a call. 1-2 times the bin size is usually appropriate.
##' @param proper Should proper mapping be counted. Default is TRUE.
##' @param map.quality the minimum mapping quality for the counted reads. Default is 30.
##' @param bp.res the base-pair resolution of the estimation. Default is 5bp.
##' @param slidW.size the size of the sliding window. Default is 200.
##' @param slidW.step the step used when sliding the window. Default is 1.
##' @return the original GRanges object with updated breakpoint location.
##' @author Jean Monlong
##' @export
breakpoint.finder <- function(bkpt.gr, files.df, test.sample, ref.samples, flank.size = c(250,2500), proper = TRUE, map.quality = 30, bp.res = c(6,8,10), slidW.size = 200, slidW.step = 1, midread=FALSE) {
 
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
    bam.df = do.call(rbind, lapply(names(bam),function(bin) data.frame(bin=bin, as.data.frame(bam[[bin]]))))
    bam.df = bam.df[which(bam.df$mapq >= map.quality),]
    return(with(bam.df, GenomicRanges::GRanges(rname, IRanges::IRanges(pos, width = qwidth), bin = bin)))
  }
  if(midread){
    cov.reads <- function(reads.gr, win.gr){
      GenomicRanges::start(reads.gr) = GenomicRanges::start(reads.gr)+GenomicRanges::width(reads.gr)
      GenomicRanges::end(reads.gr) = GenomicRanges::end(reads.gr)-GenomicRanges::width(reads.gr)+1
      GenomicRanges::countOverlaps(win.gr, reads.gr)
    }
  } else {
    cov.reads <- function(reads.gr, win.gr){
      GenomicRanges::countOverlaps(win.gr, reads.gr)
    }
  }
  comp.diff <- function(cov.m) {
    cov.m = apply(cov.m, 2, function(ee) ee/mean(ee, na.rm = TRUE))
    cov.d = median(colSums(abs(cov.m[, -1] - cov.m[, 1])), na.rm = TRUE)
    cov.d.ref = median(colSums(abs(cov.m[, -1] - cov.m[, c(3:ncol(cov.m), 2)])), na.rm = TRUE)
    return(cov.d/cov.d.ref)
  }
  middle <- function(gr) round((GenomicRanges::start(gr)+GenomicRanges::end(gr))/2)
  find.bp <- function(gr.f) {
    diff.l = lapply(bp.res, function(bpr){
      gr.bk = seq(start(gr.f), end(gr.f), bpr)
      gr.frag = GenomicRanges::GRanges(GenomicRanges::seqnames(gr.f)[1], IRanges::IRanges(gr.bk[-length(gr.bk)], width = bpr))
      cov.l = lapply(c(test.sample, ref.samples), function(samp.i) {
        gr.i = gr.bam(files.df$bam[which(files.df$sample == samp.i)], gr.f, proper = proper, map.quality = map.quality)
        cov.reads(gr.i, gr.frag)
      })
      cov = matrix(0, length(gr.frag), length(ref.samples) + 1)
      for (ii in 1:(length(ref.samples) + 1)) cov[, ii] = cov.l[[ii]]
      diff.v = sapply(seq(1, length(gr.frag) - slidW.size + 1, slidW.step), function(sw.ii) {
        comp.diff(cov[sw.ii:(sw.ii + slidW.size - 1), ])
      })
      if (all(is.na(diff.v))) {
        return(middle(gr.frag[length(gr.frag)/2]))
      } else {
        plot(diff.v)
        return(middle(gr.frag[which.max(diff.v) * slidW.step + slidW.size/2]))
      }
    })
    GenomicRanges::start(gr.f) = median(unlist(diff.l), na.rm=TRUE)
    GenomicRanges::width(gr.f) = 1
    gr.f    
  }
  
  bkpt.gr.n = lapply(bkpt.gr, function(gr.f) {
    print(gr.f)
    if(GenomicRanges::width(gr.f) > 3*flank.size[1]){
      flank.size = flank.size[2]
    }
    up.bp = find.bp(GenomicRanges::flank(gr.f, flank.size[1], both = TRUE))
    dw.bp = find.bp(GenomicRanges::flank(gr.f, flank.size[1], start = FALSE, both = TRUE))
    GenomicRanges::start(gr.f) = GenomicRanges::start(up.bp)
    GenomicRanges::end(gr.f) = GenomicRanges::end(dw.bp)
    gr.f
  })

  return(unlist(GenomicRanges::GRangesList(bkpt.gr.n)))
} 
