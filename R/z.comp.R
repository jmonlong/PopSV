##' Z-score computation from bin count and/or mean/sd metrics on the reference samples.
##'
##' The Z-score is computed by substracting the bin count by the average bin count
##' across the reference samples and dividing by their standard deviation. If
##' 'z.poisson' is TRUE, a score using Poisson distribution is also computed, using
##' the average bin count as an estimator of the lambda. Then the score with the lowest
##' absolute value is kept. This hybrid Z-score is to be used when some regions have low
##' coverage where it is more robust to use Poisson assumptions.
##' @title Z-score computation
##' @param bc.f the path to the normalized bin count file. 
##' @param norm.stats.f the name of the file with the statistic of the targeted normalization run.
##' @param files.df a data.frame with the file paths.
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @param nb.cores the number of cores to use.
##' @param chunk.size the chunk size. Default is 1e4. If NULL, no chunks are used.
##' @param append should the Z-scores be appended to existing files. Default is FALSE.
##' @return a list with
##' \item{ref.samples}{a vector with the reference samples used.}
##' \item{z.poisson}{was Normal-Poisson hybrid Z-score score computed.}
##' @author Jean Monlong
##' @export
z.comp <- function(bc.f, norm.stats.f, files.df, z.poisson = FALSE, nb.cores = 1, chunk.size=1e4, append=FALSE) {

  if(!file.exists(bc.f)){
    stop("Bin count file not found.")
  }

  if (z.poisson) {
    z.comp.f <- function(x, mean.c, sd.c) {
      z.n = (x - mean.c)/sd.c
      z.p = stats::qnorm(stats::ppois(x, mean.c))
      n.ii = abs(z.n) < abs(z.p)
      z.p[which(n.ii)] = z.n[which(n.ii)]
      z.p
    }
  } else {
    z.comp.f <- function(x, mean.c, sd.c) {
      (x - mean.c)/sd.c
    }
  }

  read.chunk <- function(chunk.start=NULL, chunk.end=NULL){
    col.n = utils::read.table(bc.f, nrows=1, sep="\t", header=FALSE, as.is=TRUE)
    bc.dt = suppressWarnings(data.table::fread(bc.f,nrows=chunk.end-chunk.start+1, skip=chunk.start, header=FALSE, sep="\t"))
    data.table::setnames(bc.dt, as.character(col.n))
    col.n = utils::read.table(norm.stats.f, nrows=1, sep="\t", header=FALSE, as.is=TRUE)
    ns.dt = suppressWarnings(data.table::fread(norm.stats.f,nrows=chunk.end-chunk.start+1, skip=chunk.start, header=FALSE, sep="\t", select=1:5))
    data.table::setnames(ns.dt, as.character(col.n)[1:5])
    list(bc=bc.dt, ns=ns.dt)
  }
  headers = as.character(utils::read.table(bc.f, nrows=1, sep="\t", header=FALSE, as.is=TRUE))
  
  ## Sample names and row number
  ref.samples = setdiff(headers, c("chr","start","end"))
  bc.1 = data.table::fread(bc.f, header = TRUE, select=1)
  nrows = nrow(bc.1)
  rm(bc.1)

  ## Compute chunk index
  if(!is.null(chunk.size) && chunk.size<nrows){
    chunks = tapply(1:nrows, rep(1:ceiling(nrows/chunk.size), each=chunk.size)[1:nrows], identity)
  } else {
    chunks = list(1:nrows)
  }

  ## For each chunk
  for(ch.ii in 1:length(chunks)){

    ## Read chunk
    chunk.o = read.chunk(min(chunks[[ch.ii]]),max(chunks[[ch.ii]]))
    bc.l = chunk.o$bc
    bc.1 = bc.l[, 1:3, with = FALSE]
    bc.l = as.matrix(bc.l[, ref.samples, with = FALSE])

    ## Get or compute mean/sd in reference
    msd = as.data.frame(chunk.o$ns)

    z = parallel::mclapply(1:ncol(bc.l), function(cc) z.comp.f(bc.l[,cc], mean.c = msd$m, sd.c = msd$sd), mc.cores=nb.cores)
    z = matrix(unlist(z), ncol=length(z))
    fc = bc.l/msd$m
    colnames(z) = colnames(fc) = ref.samples
    z = data.frame(bc.1[, 1:3, with = FALSE], z)
    fc = data.frame(bc.1[, 1:3, with = FALSE], fc)

    ## Write output files
    write.split.samples(list(z=z, fc=fc), files.df, ref.samples, files.col=c("z","fc"), compress.index=FALSE, append=append | ch.ii>1)
  }

  files.df = files.df[which(files.df$sample %in% ref.samples),]
  comp.index.files(c(files.df$z, files.df$fc), rm.input=TRUE, reorder=TRUE)

  return(list(ref.samples=ref.samples, z.poisson = z.poisson))
}
