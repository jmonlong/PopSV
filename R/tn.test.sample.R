##' Bin counts from one sample are normalized following instructions from
##' a previous targeted-normalization run.
##'
##' The Z-score is computed by substracting the bin count by the average bin count
##' across the reference samples and dividing by their standard deviation. If
##' 'z.poisson' is TRUE, a score using Poisson distribution is also computed, using
##' the average bin count as an estimator of the lambda. Then the score with the lowest
##' absolute value is kept. This hybrid Z-score is to be used when some regions have low
##' coverage where it is more robust to use Poisson assumptions.
##' @title Single sample targeted normalization and test
##' @param test.sample the name of the sample to test.
##' @param cont.sample the name of the sample used as control for the normalization.
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bc.gc.bg' are required and should be
##' present  after running 'initFileNames' function. Files should exist if
##' 'correct.GC' was run.
##' @param norm.stat.f the name of the file with the statistic of the targeted normalization run.
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @param aberrant.cases if TRUE (default) a more robust (but sligthly longer) normalization
##' is performed on cases to deal with potential large chromosomal aberrations. In practice,
##' it is recommended for cancer but can be turned off if less than ~20% of the genome is expected
##' to be affected.
##' @return a data.frame with columns :
##' \item{chr, start, end}{the location of the bin}
##' \item{bc}{the normalized bin count}
##' \item{z}{the Z-scores}
##' \item{fc}{the fold-change compared to the average bin count in the reference samples}
##' @author Jean Monlong
##' @export
tn.test.sample <- function(test.sample, cont.sample, files.df, norm.stat.f, z.poisson=FALSE, aberrant.cases=FALSE){
    chunk.size=1e3
    test.bc = read.table(subset(files.df, sample==test.sample)$bc.gc.gz, colClasses=c("character","integer","integer","numeric"), header=TRUE)
    id.test = 1:nrow(test.bc)
    names(id.test) = paste(test.bc$chr, as.integer(test.bc$start), sep="-")
    cont.bc = read.table(subset(files.df, sample==cont.sample)$bc.gc.gz, colClasses=c("character","integer","integer","numeric"), header=TRUE)
    id.cont = 1:nrow(cont.bc)
    names(id.cont) = paste(cont.bc$chr, as.integer(cont.bc$start), sep="-")
    
    if(z.poisson){
        z.comp <- function(x, mean.c, sd.c){
            z.n = (x-mean.c)/sd.c
            z.p = qnorm(ppois(x, mean.c))
            n.ii = abs(z.n) < abs(z.p)
            z.p[n.ii] = z.n[n.ii]
            z.p
        }
    } else {
        z.comp <- function(x, mean.c, sd.c){(x-mean.c)/sd.c}
    }

    res.df = createEmptyDF(c("character", rep("integer",2),rep("numeric",3)), nrow(test.bc))
    colnames(res.df) = c("chr","start","end","bc","z","fc")
    res.df$chr = test.bc$chr
    res.df$start = test.bc$start
    res.df$end = test.bc$end

    test.bin <- function(ns){
      bin = paste(ns[1],as.integer(ns[2]), sep="-")
      if(aberrant.cases){
        norm.coeff = norm.tm.opt(test.bc[id.test[ns[-(1:7)]],"bc",drop=FALSE],ref.col=cont.bc[id.cont[ns[-(1:7)]],"bc"],bc.mean.norm=as.numeric(ns[4]),chrs=test.bc[id.test[ns[-(1:7)]],"chr"])
      } else {
        norm.coeff = norm.tm.opt(test.bc[id.test[ns[-(1:7)]],"bc",drop=FALSE],ref.col=cont.bc[id.cont[ns[-(1:7)]],"bc"])
      }
      bc.n = test.bc[id.test[bin],"bc"] * norm.coeff
      return(c(bc = bc.n,
               z = z.comp(bc.n,as.numeric(ns[5]),as.numeric(ns[6])),
               fc = bc.n/as.numeric(ns[5])))
    }
    
    con = file(norm.stat.f,"r")
    headers = unlist(strsplit(readLines(con,n=1),"\t"))
    while(length((lines = readLines(con,n=chunk.size)))>0){
        norm.chunk = matrix(unlist(strsplit(lines, "\t")), ncol=length(headers), byrow=TRUE)
        colnames(norm.chunk) = headers
        norm.chunk = norm.chunk[norm.chunk[,4]!="NA",]
        bins = paste(norm.chunk[,1],as.integer(norm.chunk[,2]), sep="-")
        res.df[id.test[bins], 4:6] = t(apply(norm.chunk, 1, test.bin))
    }
    close(con)

    return(res.df)
}
