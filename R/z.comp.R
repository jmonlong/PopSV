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
##' @param recomp.msd Should the mean/SD be recomputed and written to 'norm.stats.f' ? Default is FALSE. Most of the times you don't want that but it's useful when combining discordant mapping.
##' @return a list with
##' \item{z}{a data.frame with the Z-scores if 'bc.f' was a data.frame. NULL otherwise.}
##' \item{fc}{a data.frame with the fold-changes if 'bc.f' was a data.frame. NULL otherwise.}
##' \item{msd}{a data.frame with the mean-sd information if 'bc.f' was a data.frame. NULL otherwise.}
##' \item{z.poisson}{was Normal-Poisson hybrid Z-score score computed.}
##' @author Jean Monlong
##' @export
z.comp <- function(bc.f, norm.stats.f, files.df, z.poisson = FALSE, nb.cores = 1, chunk.size=1e4, append=FALSE, recomp.msd=FALSE) {

  if (z.poisson) {
    z.comp.f <- function(x, mean.c, sd.c) {
      z.n = ifelse(sd.c<2, Inf, (x - mean.c)/sd.c)
      pv.p = stats::ppois(x, mean.c)
      z.p = stats::qnorm(ifelse(pv.p<.95, runif(length(x),0,.99), pv.p))
      n.ii = abs(z.n) < abs(z.p)
      z.p[which(n.ii)] = z.n[which(n.ii)]
      z.p
    }
    ## z.comp.f <- function(x, mean.c, sd.c) {
    ##     z.n = (x - mean.c)/sd.c
    ##     z.p = stats::qnorm(stats::ppois(x, mean.c))
    ##     n.ii = abs(z.n) < abs(z.p)
    ##     z.p[which(n.ii)] = z.n[which(n.ii)]
    ##     z.p
    ## }
  } else {
    z.comp.f <- function(x, mean.c, sd.c) {
      (x - mean.c)/sd.c
    }
  }

  if(is.data.frame(bc.f)){
    read.chunk <- function(){
      if(!firstChunk) return(NULL)
      list(bc=bc.f, ns=norm.stats.f)
    }
  } else {
    if(!file.exists(bc.f)){
      stop("Bin count file not found.")
    }

    ## Open file connection and prepare which columns to read
    con.bc = file(bc.f, "r")
    bc.header = unlist(strsplit(readLines(con.bc, n = 1), "\t"))
    bc.colClasses = c("character", rep("numeric",2), ifelse(bc.header[-(1:3)] %in% files.df$sample, "numeric", "NULL"))
    bc.header = bc.header[which(bc.colClasses != "NULL")]
    files.df = files.df[which(files.df$sample %in% bc.header),]
    if(!recomp.msd){
      con.ns = file(norm.stats.f, "r")
      ns.header = unlist(strsplit(readLines(con.ns, n = 1), "\t"))
      ns.colClasses = c("character", rep("numeric",2), ifelse(ns.header[-(1:3)] %in% c("m","sd"), "numeric", "NULL"))
      ns.header = ns.header[1:5]
    }

    read.chunk <- function(){
      bc.res = tryCatch(utils::read.table(con.bc, colClasses=bc.colClasses, nrows=chunk.size), error=function(e)return(NULL))
      if(is.null(bc.res)){
        return(NULL)
      }
      colnames(bc.res) = bc.header
      if(!recomp.msd){
        ns.res = tryCatch(utils::read.table(con.ns, colClasses=ns.colClasses, nrows=chunk.size), error=function(e)return(NULL))
        if(is.null(ns.res)){
          stop("Different number of rows between ", bc.f," and ", norm.stats.f, " !?")
        }
        colnames(ns.res) = ns.header
      } else {
        ns.res = NULL
      }
      list(bc=bc.res, ns=ns.res)
    }
  }

  ## For each chunk until it's NULL
  firstChunk = TRUE
  while(!is.null((chunk.o = read.chunk()))){

    ## Read chunk
    bc.l = chunk.o$bc
    bc.1 = bc.l[, 1:3]
    bc.l = as.matrix(bc.l[, -(1:3)])

    if(recomp.msd){
      msd = parallel::mclapply(1:nrow(bc.l), function(rr) unlist(mean.sd.outlierR(bc.l[rr,])), mc.cores=nb.cores)
      msd = matrix(unlist(msd), nrow=3)
      rownames(msd) = c("m","sd","nb.remove")
      msd = cbind(bc.1, t(msd))
      if(!is.data.frame(bc.f)){
        utils::write.table(msd, file=norm.stats.f, row.names=FALSE, sep="\t", col.names=firstChunk, append=!firstChunk)
      }
    } else {
      msd = chunk.o$ns
    }

    z = parallel::mclapply(1:ncol(bc.l), function(cc) z.comp.f(bc.l[,cc], mean.c = msd$m, sd.c = msd$sd), mc.cores=nb.cores)
    z = matrix(unlist(z), ncol=length(z))
    fc = bc.l/msd$m
    colnames(z) = colnames(fc) = colnames(bc.l)
    z = data.frame(bc.1, z, stringsAsFactors=TRUE)
    fc = data.frame(bc.1, fc, stringsAsFactors=TRUE)

    ## Write output files
    if(!is.data.frame(bc.f)){
      write.split.samples(list(z=z, fc=fc), files.df, files.col=c("z","fc"), compress.index=FALSE, append=append | !firstChunk)
      z = fc = msd = NULL
    }
    firstChunk = FALSE
  }

  if(!is.data.frame(bc.f)){
    close(con.bc)
    if(!recomp.msd){
      close(con.ns)
    }
    comp.index.files(c(files.df$z, files.df$fc), rm.input=TRUE, reorder=TRUE)
  }

  return(list(z=z, fc=fc, msd=msd, z.poisson = z.poisson))
}
