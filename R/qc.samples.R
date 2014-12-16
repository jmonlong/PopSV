##' Join bin counts of the samples to analyze and compute some QC metrics to help
##' defining the set of reference samples/
##' @title Join and QC of samples
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bc.gc.bg' are required and should be
##' present  after running 'initFileNames' function. Files should exist if
##' 'correct.GC' was run.
##' @param bin.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param ref.samples a vector with the names of the samples planned to be used
##' as reference.
##' @param outfile.prefix the prefix of the name of the output file. The suffix '.bgz' will
##' be appended to this name prefix after compression. 
##' @param out.pdf the name of the output pdf file.
##' @param appendIndex.outfile if TRUE (default), the results will be appended regularly on
##' the output file which will be ultimately indexed. This is recommend when a large number
##' of bins are analyzed. If FALSE, a data.frame with the bin counts will be returned and no
##' file are created.
##' @param chunk.size the number of bins to analyze at a time (for memory optimization).
##' Default is 100 000. Reduce this number if memory problems arise.
##' @param col.bc the column from 'files.df' defining the bin count file names.
##' @param nb.cores number of cores to use. If higher than 1, \code{parallel}
##' package is used to parallelized the counting.
##' @param nb.bin.support the number of supporting bins to define for the future normalization.
##' @return a list with
##' \item{bc}{the name of the file with the joined bin counts OR a data.frame with
##' these bin counts.}
##' \item{dstat}{a data.frame with for each sample its D-statistic (average correlation
##' with other reference samples).}
##' \item{pc.1.3}{a matrix with the first 3 principal components for all samples.}
##' \item{cor.pw}{a matrix with correlation for any pair of samples.}
##' @author Jean Monlong
##' @export
qc.samples <- function(files.df, bin.df, ref.samples=NULL, outfile.prefix, out.pdf=NULL, appendIndex.outfile=TRUE, chunk.size=1e5, col.bc="bc.gc.gz", nb.cores=1, nb.bin.support=1e4){
    if(is.null(ref.samples)) ref.samples = as.character(files.df$sample)
    if(nrow(bin.df)<1.3*chunk.size){
        bc.df = createEmptyDF(c("character",rep("integer",2), rep("numeric",nrow(files.df))), nrow(bin.df))
        colnames(bc.df) = c("chr","start","end", as.character(files.df$sample))
        bc.df$chr = bin.df$chr
        bc.df$start = bin.df$start
        bc.df$end = bin.df$end
        if(nb.cores>1){
            bc.l = parallel::mclapply(files.df[,col.bc], function(fi){
                read.bedix(fi, bin.df)[,4]
            },mc.cores=nb.cores)
        } else {
            bc.l = lapply(files.df[,col.bc], function(fi){
                read.bedix(fi, bin.df)[,4]
            })
        }
        med.samp = lapply(bc.l,median, na.rm=TRUE) ## Normalization of the median coverage
        med.med = median(unlist(med.samp), na.rm=TRUE)
        for(samp.i in 1:nrow(files.df)){
            bc.df[,as.character(files.df$sample[samp.i])] = bc.l[[samp.i]] * med.med / med.samp[[samp.i]]
        }
        med.cov.df = data.frame(bc.df[,1:3], med.bc= apply(bc.df[,ref.samples], 1, median, na.rm=TRUE))
        write.table(bc.df, file=outfile.prefix, quote=FALSE, row.names=FALSE, sep="\t")
    } else {
        nb.chunks = ceiling(nrow(bin.df)/chunk.size)
        bin.df$chunk = rep(1:nb.chunks,each=chunk.size)[1:nrow(bin.df)]
        ## Median coverage
        med.samp = rep(list(1),nrow(files.df)) ## Normalization of the median coverage
        med.med = 1
        analyze.chunk <- function(df, write.out=TRUE, sub.bc=TRUE){
          ch.nb = as.numeric(df$chunk[1])
          bc.df = createEmptyDF(c("character",rep("integer",2), rep("numeric",nrow(files.df))), nrow(df))
          colnames(bc.df) = c("chr","start","end", as.character(files.df$sample))
          bc.df$chr = df$chr
          bc.df$start = df$start
          bc.df$end = df$end
          if(nb.cores>1){
            bc.l = parallel::mclapply(files.df[,col.bc], function(fi){
              read.bedix(fi, df)[,4]
            },mc.cores=nb.cores)
          } else {
            bc.l = lapply(files.df[,col.bc], function(fi){
              read.bedix(fi, df)[,4]
            })
          }
          for(samp.i in 1:nrow(files.df)){
            bc.df[,as.character(files.df$sample[samp.i])] = bc.l[[samp.i]] * med.med / med.samp[[samp.i]]
          }
          med.cov.df = data.frame(bc.df[,1:3], med.bc= apply(bc.df[,ref.samples], 1, median, na.rm=TRUE))
          if(write.out){
            write.table(bc.df, file=outfile.prefix, quote=FALSE, row.names=FALSE, sep="\t", append=ch.nb>1, col.names=ch.nb==1)
          }
          if(sub.bc) bc.df=bc.df[sample(1:nrow(bc.df),chunk.size/nb.chunks),]
          return(list(med.cov.df=med.cov.df, bc.df=bc.df))
        }
        bc.rand = analyze.chunk(bin.df[sample(1:nrow(bin.df), 1e3),], write.out=FALSE, sub.bc=FALSE)
        med.samp = lapply(as.character(files.df$sample), function(samp.i)median(bc.rand$bc.df[,samp.i], na.rm=TRUE)) ## Normalization of the median coverage
        med.med = median(unlist(med.samp), na.rm=TRUE)
        bc.res = lapply(unique(bin.df$chunk), function(chunk.i){
            analyze.chunk(subset(bin.df, chunk==chunk.i))
        })
        med.cov.df = plyr::ldply(bc.res, function(ee)ee$med.cov.df)
        bc.df = plyr::ldply(bc.res, function(ee)ee$bc.df)
    }

    ## Select set of supporting bins
    med.cov.df = subset(med.cov.df, med.bc!=0)
    q.bk = quantile(med.cov.df$med.bc, c(0,.05,.1,.2,.8,.9,.95,1))
    med.cov.cut = cut(med.cov.df$med.bc, breaks=q.bk)
    cut.prop = c(.2,.15,.1,.1,.1,.15,.2)/2
    bin.sup.i = unique(unlist(sapply(1:nlevels(med.cov.cut), function(ii)sample(which(as.numeric(med.cov.cut)==ii), nb.bin.support*cut.prop[ii], replace=TRUE))))
    bin.sup.i = c(bin.sup.i, sample(setdiff(1:nrow(med.cov.df), bin.sup.i), nb.bin.support-length(bin.sup.i)))
    bin.sup.df = med.cov.df[bin.sup.i, ]
    ##

    ## Compress and index
    if(appendIndex.outfile){
        final.file = paste(outfile.prefix,".bgz",sep="")
        Rsamtools::bgzip(outfile.prefix, dest=final.file, overwrite=TRUE)
        file.remove(outfile.prefix)
        Rsamtools::indexTabix(final.file, format="bed")
    }
    ##
    
    ## QC
    cor.bs <- function(x,nbsim=10,prop=.1,replace=FALSE){
        ni.sim = nrow(x)*prop
        cor.sim = sapply(1:nbsim,function(i)as.vector(cor(x[sample(1:nrow(x),ni.sim,replace=replace),],use="comp")))
        res = matrix(apply(cor.sim,1,median),ncol(x),ncol(x))
        colnames(res) = rownames(res) = colnames(x)
        return(res)
    }
    all.samples = colnames(bc.df)[-(1:3)]
    bc.mv = medvar.norm.internal(bc.df[,all.samples])
    corbs = cor.bs(bc.mv)
    cor.pw = corbs
    diag(cor.pw) = NA
    meanCor = apply(cor.pw[ref.samples,ref.samples],1,function(e)mean(e,na.rm=TRUE))
    ## PCA
    pc = prcomp(t(na.exclude(bc.mv)))
    if(!is.null(out.pdf)){
        pdf(out.pdf,10,8)
        hc = hclust(as.dist(1-corbs))
        plot(hc)
        print(ggplot2::qplot(x=meanCor) +
              ggplot2::geom_histogram() +
              ggplot2::ggtitle("D statistic - Mean correlation with the other samples") +
              ggplot2::theme_bw())
        plot(pc$x[,1:2],type="n",main="PCA")
        text(pc$x[,1:2],labels=all.samples)
        dev.off()
    }
    if(appendIndex.outfile){
        bc.df = paste(outfile.prefix,".bgz",sep="")
    }
    return(list(bc=bc.df,dstat=data.frame(sample=ref.samples,Dstat=meanCor),
                pc.1.3=pc$x[,1:3], cor.pw=corbs, bin.sup.df=bin.sup.df))
}
