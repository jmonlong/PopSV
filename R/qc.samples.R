##' Join bin counts of the samples to analyze and compute some QC metrics to help
##' defining the set of reference samples.
##' @title Join and QC the reference samples
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bc.gc.bg' are required and should be
##' present  after running 'initFileNames' function. Files should exist if
##' 'correct.GC' was run.
##' @param bin.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param outfile.prefix the prefix of the output file name. The suffix '.bgz' will
##' be appended if compressed ('appendIndex.outfile=TRUE').
##' @param ref.samples a vector with the names of the samples to use as reference.
##' @param nb.ref.samples the number of reference samples desired. If NULL, the size of 'ref.samples'.
##' @param plot should PCA graphs be outputed ? Default is TRUE.
##' @param appendIndex.outfile if TRUE (default), the results will be appended regularly on
##' the output file which will be ultimately compressed and indexed. This is recommend when a large number
##' of bins are analyzed. If FALSE, a data.frame with the bin counts will be returned and no
##' file are created.
##' @param chunk.size the number of bins to analyze at a time (for memory optimization).
##' Default is 100 000. Reduce this number if memory problems arise.
##' @param col.bc the column from 'files.df' defining the bin count file names.
##' @param nb.cores number of cores to use. If higher than 1, \code{parallel}
##' package is used to parallelize the counting.
##' @return a list with
##' \item{bc}{the name of the file with the joined bin counts OR a data.frame with
##' these bin counts.}
##' \item{ref.samples}{a vector with the reference samples names.}
##' \item{cont.sample}{the name of the sample to use as control among the reference samples (for normalization).}
##' \item{pc.all.df}{a data.frame with the first 3 principal components for all input reference samples.}
##' \item{pc.ref.df}{a data.frame with the first 3 principal components for the final reference samples.}
##' @author Jean Monlong
##' @export
qc.samples <- function(files.df, bin.df, outfile.prefix, ref.samples = NULL, nb.ref.samples = NULL, plot = TRUE, appendIndex.outfile = TRUE, chunk.size = 1e+05, col.bc = "bc.gc.gz", nb.cores = 1) {

  ## Checks and arguments completion
  if(!all(c("chr","start","end") %in% colnames(bin.df))){
    stop("Missing column in 'bin.df'. 'chr', 'start' and 'end' are required.")
  }
  if(!all(c("sample", col.bc) %in% colnames(files.df))){
    stop("Missing column in 'files.df'. 'sample', '",col.bc,"' are required.")
  }
  if (is.null(ref.samples)) {
    ref.samples = as.character(files.df$sample)
  }
  if (is.null(nb.ref.samples)) {
    nb.ref.samples = length(ref.samples)
  }

  ## Flexible function to read bin counts on specific bins, several files, ...
  read.bc.samples <- function(df, files.df, nb.cores = 1, med.med = NULL, med.samp = NULL,
                              file = NULL, append.f = FALSE, sub.bc = NULL) {
    bc.df = createEmptyDF(c("character", rep("integer", 2), rep("numeric", nrow(files.df))),
      nrow(df))
    colnames(bc.df) = c("chr", "start", "end", as.character(files.df$sample))
    bc.df$chr = df$chr
    bc.df$start = df$start
    bc.df$end = df$end
    if (nb.cores > 1) {
      bc.l = parallel::mclapply(as.character(files.df[, col.bc]), function(fi) {
        read.bedix(fi, df,exact.match=TRUE)[, 4]
      }, mc.cores = nb.cores)
    } else {
      bc.l = lapply(as.character(files.df[, col.bc]), function(fi) {
        read.bedix(fi, df,exact.match=TRUE)[, 4]
      })
    }
    if (is.null(med.med) & is.null(med.samp)) {
      med.samp = lapply(bc.l, median, na.rm = TRUE)  ## Normalization of the median coverage
      med.med = median(unlist(med.samp), na.rm = TRUE)
    }
    for (samp.i in 1:nrow(files.df)) {
      bc.df[, as.character(files.df$sample[samp.i])] = round(bc.l[[samp.i]] * med.med/med.samp[[samp.i]], 2)
    }
    bc.df = bc.df[order(bc.df$chr, bc.df$start, bc.df$end),]
    if (!is.null(file)) {
      write.table(bc.df, file = file, quote = FALSE, row.names = FALSE, sep = "\t",
                  append = append.f, col.names = !append.f)
    }
    if (!is.null(sub.bc)) {
      bc.df = bc.df[sample(1:nrow(bc.df), min(c(nrow(bc.df), sub.bc))), ]
    }
    return(bc.df)
  }
  center.pt <- function(m) {
    dm = as.matrix(dist(m))
    diag(dm) = NA
    rownames(m)[which.min(apply(dm, 2, mean, na.rm = TRUE))]
  }

  files.df = files.df[which(files.df$sample %in% ref.samples), ]
  ## Too many ref.samples
  pc.all.df = bc.rand = NULL
  if (length(ref.samples) > nb.ref.samples) {
    sparse.pts <- function(m, nb.pts) {
      dm = as.matrix(dist(m))
      diag(dm) = NA
      pts.jj = 1:nrow(dm)
      pts.ii = which.max(apply(dm, 1, max, na.rm = TRUE))
      pts.jj = pts.jj[-pts.ii]
      while (length(pts.ii) < nb.pts) {
        m.ii = which.max(apply(dm[pts.ii, pts.jj, drop = FALSE], 2, min,
          na.rm = TRUE))
        pts.ii = c(pts.ii, pts.jj[m.ii])
        pts.jj = pts.jj[-m.ii]
      }
      m[pts.ii, ]
    }
    ## Get subset of bins
    bc.rand = read.bc.samples(bin.df[sample.int(nrow(bin.df), min(c(nrow(bin.df)/2,1000))), ], files.df, med.med = 1, med.samp = rep(list(1), nrow(files.df)))
    ## PCA
    bc.mv = medvar.norm.internal(bc.rand[, ref.samples])
    pc.all = prcomp(t(na.exclude(bc.mv)))
    pc.all.df = data.frame(pc.all$x[, 1:3])
    pc.all.df$sample = ref.samples
    sp.o = sparse.pts(pc.all$x[, 1:2], nb.ref.samples)
    pc.all.df$reference = pc.all.df$sample %in% rownames(sp.o)
    if (plot) {
      PC1 = PC2 = reference = NULL  ## Uglily appease R checks
      print(ggplot2::ggplot(pc.all.df, ggplot2::aes(x = PC1, y = PC2, size = reference)) +
              ggplot2::geom_point(alpha = 0.8) + ggplot2::theme_bw() + ggplot2::scale_size_manual(values = 2:3))
    }
    ref.samples = rownames(sp.o)
    bc.rand = bc.rand[, c("chr", "start", "end", ref.samples)]
  }

  files.df = files.df[which(files.df$sample %in% ref.samples), ]
  bin.df = bin.df[order(bin.df$chr, bin.df$start, bin.df$end),]
  if (nrow(bin.df) < 1.3 * chunk.size) {
    bc.df = read.bc.samples(bin.df, files.df, nb.cores)
    write.table(bc.df, file = outfile.prefix, quote = FALSE, row.names = FALSE,
                sep = "\t")
  } else {
    nb.chunks = ceiling(nrow(bin.df)/chunk.size)
    bin.df$chunk = rep(1:nb.chunks, each = chunk.size)[1:nrow(bin.df)]
    analyze.chunk <- function(df) {
      ch.nb = as.numeric(df$chunk[1])
      df.o = read.bedix(as.character(files.df[1, col.bc]), df,exact.match=TRUE)
      read.bc.samples(df.o, files.df, nb.cores, med.med, med.samp, file = outfile.prefix,
                      append.f = ch.nb > 1, sub.bc = chunk.size/nb.chunks)
    }
    if (is.null(bc.rand)) {
      bc.rand = read.bc.samples(bin.df[sample.int(nrow(bin.df), min(c(nrow(bin.df)/2,1000))), ], files.df, med.med = 1, med.samp = rep(list(1), nrow(files.df)))
    }
    med.samp = lapply(as.character(files.df$sample), function(samp.i) median(bc.rand[,
      samp.i], na.rm = TRUE))  ## Normalization of the median coverage
    med.med = median(unlist(med.samp), na.rm = TRUE)
    bc.res = lapply(unique(bin.df$chunk), function(chunk.i) {
      analyze.chunk(bin.df[which(bin.df$chunk == chunk.i), ])
    })
    bc.df = plyr::ldply(bc.res, identity)
  }

  ## Compress and index
  if (appendIndex.outfile) {
    outfile.prefix = comp.index.files(outfile.prefix, rm.input=TRUE, overwrite.out=TRUE, reorder=any(order(bin.df$chr, bin.df$start, bin.df$end)!=1:nrow(bin.df)))
  }
  ##

  ## QC
  bc.mv = medvar.norm.internal(bc.df[, ref.samples])
  if(any(is.infinite(bc.mv))){
    bc.mv[which(is.infinite(bc.mv))] = NA
  }
  pc.ref = prcomp(t(na.exclude(bc.mv)))
  pc.ref.df = data.frame(pc.ref$x[, 1:3])
  pc.ref.df$sample = ref.samples
  if (plot) {
    PC1 = PC2 = NULL  ## Uglily appease R checks
    print(ggplot2::ggplot(pc.ref.df, ggplot2::aes(x = PC1, y = PC2)) + ggplot2::geom_point(alpha = 0.8) +
            ggplot2::theme_bw())
  }
  cont.sample = center.pt(pc.ref$x[, 1:2])
  bc.df = outfile.prefix
  return(list(bc = bc.df, ref.samples = ref.samples, cont.sample = cont.sample,
              pca.all.df = pc.all.df, pca.ref.df = pc.ref.df))
}
