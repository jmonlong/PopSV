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
##' @return a list with
##' \item{z}{a data.frame with the Z-scores for each bin and sample (bin x sample).}
##' \item{fc}{a data.frame with the fold-change compared to the average bin count in
##' the reference samples for each bin and sample (bin x sample).}
##' \item{msd}{the mean, standard deviation and number of removed outlier samples in each bin.}
##' @author Jean Monlong
##' @export
z.comp <- function(bc.f, files.df, ref.samples=NULL, z.poisson = FALSE, nb.cores = 1, chunk.size=NULL, out.msd.f="ref-msd.tsv") {

  if(!file.exists(bc.f)){
    stop("Bin count file not found.")
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

  read.chunk <- function(chunk.start=NULL, chunk.end=NULL, file, sep="\t", header=TRUE){
    if(header){
      col.n = read.table(file, nrows=1, sep=sep, header=FALSE, as.is=TRUE)
    }
    dt = suppressWarnings(data.table::fread(file,nrows=chunk.end-chunk.start+1, skip=chunk.start-1+as.numeric(header), header=FALSE, sep=sep))
    data.table::setnames(dt, as.character(col.n))
    dt
  }

  ## One column to get number of rows
  bc.1 = data.table::fread(bc.f, header = TRUE, select=1)
  nrows = nrow(bc.1)
  rm(bc.1)

  ## Sample names
  headers = as.character(read.table(bc.f, nrows=1, sep="\t", header=FALSE, as.is=TRUE))
  if(!is.null(ref.samples)){
    if(!all(ref.samples %in% headers)){
      stop("Some specified reference samples are not present in the bin count file.")
    }
  } else {
    ref.samples = setdiff(headers, c("chr","start","end"))
  }
  
  ## Compute chunk index
  if(!is.null(chunk.size)){
    chunks = tapply(1:nrows, rep(1:ceiling(nrows/chunk.size), each=chunk.size)[1:nrows], identity)
  } else {
    chunks = list(1:nrows)
  }
   
  ## For each chunk
  for(ch.ii in 1:length(chunks)){

    ## Read chunk
    bc.l = read.chunk(min(chunks[[ch.ii]]),max(chunks[[ch.ii]]),bc.f)
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
    if(!is.null(chunk.size)){
      write.split.samples(list(z=z, fc=fc), files.df, ref.samples, res.n=c("z","fc"), files.col=c("z","fc"), compress.index=FALSE, append=ch.ii>1)
    }
    
    ## Write mean/sd file
    msd = data.frame(as.data.frame(bc.1[, 1:3, with=FALSE]), t(msd))
    if(!is.null(chunk.size) & !is.null(out.msd.f)){
      write.table(msd, file=out.msd.f, sep="\t", row.names=FALSE, quote=FALSE, append=ch.ii>1, col.names=ch.ii==1)
    }
    
  }
  
  return(list(z = z, fc = fc, msd = msd, z.poisson = z.poisson))
}

## To test: one sample only; chunks; msd input or not 
## To add: Check consistency in parameters in the beginning
