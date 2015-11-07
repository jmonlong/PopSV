##' Z-score computation from bin count and/or mean/sd metrics on the reference samples
##'
##' @title Z-score computation
##' @param bc.f the path to the normalized bin count file.
##' @param files.df a data.frame with the file paths.
##' @param ref.samples a vector with the samples to use as reference. If 'NULL' (default) all samples in the bin count file are used as reference.
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @param nb.cores the number of cores to use.
##' @param chunk.size the chunk size. If NULL (Default), no chunks are used.
##' @param out.msd.f the name of the file to write the mean/sd information. If 'NULL' no file is created.
##' @param append should the Z-scores be appended to existing files. Default is FALSE.
##' @param files.col the name of the column from 'files.df' to use to get the bin counts. Used only if 'bc.f' is NULL.
##' @return a list with
##' \item{ref.samples}{a vector with the reference samples used.}
##' \item{z.poisson}{was Normal-Poisson hybrid Z-score score computed.}
##' @author Jean Monlong
##' @export
z.comp <- function(bc.f=NULL, files.df, ref.samples=NULL, z.poisson = FALSE, nb.cores = 1, chunk.size=NULL, out.msd.f="ref-msd.tsv", append=FALSE, files.col="bc.gc.norm.gz") {

  if(!is.null(bc.f) && !file.exists(bc.f)){
    stop("Bin count file not found.")
  }
  if(is.null(bc.f) & is.null(files.col)){
    stop("Either 'bc.f' or 'files.col' must be non-NULL")
  }
  if(is.null(bc.f) && !is.null(files.col) && all(grepl("\\.bgz$", files.df[,files.col]))){
    files.df[,files.col] = paste("zcat",files.df[,files.col])
  }

  if (z.poisson) {
    z.comp.f <- function(x, mean.c, sd.c) {
      z.n = (x - mean.c)/sd.c
      z.p = qnorm(ppois(x, mean.c))
      n.ii = abs(z.n) < abs(z.p)
      z.p[which(n.ii)] = z.n[which(n.ii)]
      z.p
    }
  } else {
    z.comp.f <- function(x, mean.c, sd.c) {
      (x - mean.c)/sd.c
    }
  }

  if(!is.null(bc.f)){
    read.chunk <- function(chunk.start=NULL, chunk.end=NULL){
      col.n = read.table(bc.f, nrows=1, sep="\t", header=FALSE, as.is=TRUE)
      dt = suppressWarnings(data.table::fread(bc.f,nrows=chunk.end-chunk.start+1, skip=chunk.start, header=FALSE, sep="\t"))
      data.table::setnames(dt, as.character(col.n))
      dt
    }
    headers = as.character(read.table(bc.f, nrows=1, sep="\t", header=FALSE, as.is=TRUE))
  } else {
    read.chunk <- function(chunk.start=NULL, chunk.end=NULL){
      col.n = c("chr","start","end", files.df$sample)
      files = as.character(files.df[,files.col])
      dt.l = parallel::mclapply(files, function(file){
        suppressWarnings(data.table::fread(file,nrows=chunk.end-chunk.start+1, skip=chunk.start, header=FALSE, sep="\t"))
      }, mc.cores = nb.cores)
      bc.1 = dt.l[[1]][, 1:3, with = FALSE]
      if(any(!unlist(parallel::mclapply(dt.l, function(dt) all(bc.1[,2,with=FALSE]==dt[,2,with=FALSE] & bc.1[,1,with=FALSE]==dt[,1,with=FALSE]), mc.cores = nb.cores)))){
        stop("Different orders in the bin count files.")
      }
      dt = do.call(cbind, lapply(dt.l, function(dt) dt[,4,with=FALSE]))
      dt = cbind(bc.1, dt)
      data.table::setnames(dt, as.character(col.n))
      dt
    }
    headers = as.character(files.df$sample)
  }

  ## Sample names and row number
  if(!is.null(ref.samples)){
    if(!all(ref.samples %in% headers)){
      stop("Some specified reference samples are not present in the bin count file.")
    }
  } else {
    ref.samples = setdiff(headers, c("chr","start","end"))
  }
  if(!is.null(bc.f)){
    bc.1 = data.table::fread(bc.f, header = TRUE, select=1)
  } else {
    files.df = subset(files.df, sample %in% ref.samples)
    bc.1 = data.table::fread(as.character(files.df[,files.col])[1], header = TRUE, select=1)
  }
  nrows = nrow(bc.1)
  rm(bc.1)

  ## Compute chunk index
  if(!is.null(chunk.size)){
    chunks = tapply(1:nrows, rep(1:ceiling(nrows/chunk.size), each=chunk.size)[1:nrows], identity)
  } else {
    chunks = list(1:nrows)
  }

  ## For each chunk
  for(ch.ii in 1:length(chunks)){

    ## Read chunk
    bc.l = read.chunk(min(chunks[[ch.ii]]),max(chunks[[ch.ii]]))
    bc.1 = bc.l[, 1:3, with = FALSE]
    bc.l = as.matrix(bc.l[, ref.samples, with = FALSE])

    ## Get or compute mean/sd in reference
    msd = parallel::mclapply(1:nrow(bc.l), function(rr) unlist(mean.sd.outlierR(bc.l[rr,])), mc.cores=nb.cores)
    msd = matrix(unlist(msd), nrow=3)
    rownames(msd) = c("m","sd","nb.remove")

    z = parallel::mclapply(1:ncol(bc.l), function(cc) z.comp.f(bc.l[,cc], mean.c = msd[1, ], sd.c = msd[2, ]), mc.cores=nb.cores)
    z = matrix(unlist(z), ncol=length(z))
    fc = bc.l/msd[1, ]
    colnames(z) = colnames(fc) = ref.samples
    z = data.frame(bc.1[, 1:3, with = FALSE], z)
    fc = data.frame(bc.1[, 1:3, with = FALSE], fc)

    ## Write output files
    write.split.samples(list(z=z, fc=fc), files.df, ref.samples, files.col=c("z","fc"), compress.index=FALSE, append=append | ch.ii>1)

    ## Write mean/sd file
    msd = data.frame(as.data.frame(bc.1[, 1:3, with=FALSE]), t(msd))
    if(!is.null(out.msd.f)){
      write.table(msd, file=out.msd.f, sep="\t", row.names=FALSE, quote=FALSE, append=append | ch.ii>1, col.names=!append & ch.ii==1)
    }

  }

  files.df = files.df[which(files.df$sample %in% ref.samples),]
  comp.index.files(c(files.df$z, files.df$fc), rm.input=TRUE, reorder=TRUE)

  return(list(ref.samples=ref.samples, z.poisson = z.poisson))
}

## To test: one sample only; chunks; msd input or not
## To add: Check consistency in parameters in the beginning
