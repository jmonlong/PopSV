##' An attempt to find breakpoints as the location with the most dramatic local coverage differences between a samples and a set of reference samples. DON'T TRUST THIS FUNCTION YET.
##' @title Breakpoint finder (NOT READY YET)
##' @param range.df a data.frame with the information about the genomic ranges to analyze.
##' @param files.df a data.frame with information about the BAM file path for each sample (columns 'sample' and 'bam').
##' @param test.sample the name of the sample of interest.
##' @param ref.samples the names of the samples to use as reference for the breakpoint detection.
##' @param proper Should proper mapping be counted. Default is TRUE.
##' @param map.quality the minimum mapping quality for the counted reads. Default is 30.
##' @param bp.res the base-pair resolution of the estimation.
##' @param slidW.step the step used when sliding the window. Default is 1.
##' @param midread should the middle of the read used to asses coverage. If FALSE (default), the entire read counts for coverage computation.
##' @param nb.cores the number of processing cores to use, Default is 1.
##' @return the original data.frame with updated breakpoint location.
##' @author Jean Monlong
##' @export
breakpoint.finder <- function(range.df, files.df, test.sample, ref.samples, proper = TRUE, map.quality = 0, bp.res = 3:6, slidW.step = 1, midread=FALSE, nb.cores=1) {

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
    cor(cov.m[,1], cov.m[,-1])
  }
  middle <- function(gr) round((GenomicRanges::start(gr)+GenomicRanges::end(gr))/2)
  find.bp <- function(chr, start, end, nb.bin.cons, upstream=TRUE) {
    if(nb.bin.cons==1){
      min.size = (end-start+1)/2
    } else {
      min.size = (end-start+1)/nb.bin.cons
    }
    if(upstream){
      end = start + min.size
      start = start - min.size
    } else {
      start = end - min.size
      end = end + min.size
    }
    slidW.size = floor(min(200, min.size/max(bp.res)))
    gr.f = GenomicRanges::GRanges(chr, IRanges::IRanges(start-min.size/2, end+min.size/2))
    reads.l = lapply(c(test.sample, ref.samples), function(samp.i) {
      gr.bam(files.df$bam[which(files.df$sample == samp.i)], gr.f, proper = proper, map.quality = map.quality)
    })
    diff.l = lapply(bp.res, function(bpr){
      start = start - (slidW.size*bpr/2)
      end = end + (slidW.size*bpr/2)
      gr.bk = seq(start, end, bpr)
      gr.frag = GenomicRanges::GRanges(chr, IRanges::IRanges(gr.bk[-length(gr.bk)], width = bpr))
      cov.l = lapply(1:length(reads.l), function(ii) {
        cov.reads(reads.l[[ii]], gr.frag)
      })
      cov = matrix(0, length(gr.frag), length(ref.samples) + 1)
      for (ii in 1:(length(ref.samples) + 1)) cov[, ii] = cov.l[[ii]]
      diff.m = sapply(seq(1, length(gr.frag) - slidW.size + 1, slidW.step), function(sw.ii) {
        comp.diff(cov[sw.ii:(sw.ii + slidW.size - 1), ])
      })
      apply(diff.m, 1, function(diff.v){
        if (all(is.na(diff.v))) {
          return(middle(gr.frag[length(gr.frag)/2]))
        } else {
          ##plot(diff.v)
          return(middle(gr.frag[which.min(diff.v) * slidW.step + slidW.size/2]))
        }
      })
    })
    list(bp = median(unlist(diff.l), na.rm=TRUE), sd = sd(unlist(diff.l), na.rm=TRUE))
  }

  if(all(colnames(range.df)!="nb.bin.cons")){
    warning("No 'nb.bin.cons' column found. For safety single-bin call is assumed.")
    range.df$nb.bin.cons = 1
  }

  range.df$former.size = with(range.df, end-start+1)
  res.l = parallel::mclapply(1:nrow(range.df), function(ii) {
    up.bp = find.bp(range.df$chr[ii], range.df$start[ii], range.df$end[ii], range.df$nb.bin.cons[ii])
    dw.bp = find.bp(range.df$chr[ii], range.df$start[ii], range.df$end[ii], range.df$nb.bin.cons[ii], upstream=FALSE)
    list(start=up.bp$bp, end=dw.bp$bp, bpsd.up=up.bp$sd, bpsd.dw=dw.bp$sd)
  }, mc.cores=nb.cores)

  bkres.df = range.df
  bkres.df$former.start = bkres.df$start
  bkres.df$former.end = bkres.df$end
  bkres.df$start = unlist(lapply(res.l, function(ee)ee$start))
  bkres.df$end = unlist(lapply(res.l, function(ee)ee$end))
  bkres.df$bpsd.up = unlist(lapply(res.l, function(ee)ee$bpsd.up))
  bkres.df$bpsd.dw = unlist(lapply(res.l, function(ee)ee$bpsd.dw))
  bkres.df$bp.in = with(bkres.df, apply(cbind(former.end-start,end-former.start,former.end-former.start,end-start),1,min))+1
  bkres.df$size = with(bkres.df, end-start+1)
  bkres.df$bpsd.m = with(bkres.df, (bpsd.up+bpsd.dw)/2)
  bkres.df$cn = 2*bkres.df$fc
  bkres.df$cn.new = with(bkres.df, (former.size*cn-2*(former.size-bp.in))/bp.in)
  return(bkres.df)
}
