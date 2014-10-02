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
##' @return a list with
##' \item{bc}{the name of the file with the joined bin counts OR a data.frame with
##' these bin counts.}
##' \item{dstat}{a data.frame with for each sample its D-statistic (average correlation
##' with other reference samples).}
##' \item{pc.1.3}{a matrix with the first 3 principal components for all samples.}
##' \item{cor.pw}{a matrix with correlation for any pair of samples.}
##' @author Jean Monlong
##' @export
qc.samples <- function(files.df, bin.df, ref.samples=NULL, outfile.prefix, out.pdf=NULL, appendIndex.outfile=TRUE, chunk.size=1e5){
    if(nrow(bin.df)<1.3*chunk.size){
        bc.df = createEmptyDF(c("character",rep("integer",2), rep("numeric",nrow(files.df))), nrow(bin.df))
        colnames(bc.df) = c("chr","start","end", as.character(files.df$sample))
        bc.df$chr = bin.df$chr
        bc.df$start = bin.df$start
        bc.df$end = bin.df$end
        for(samp.i in 1:nrow(files.df)){
            bc.df[,as.character(files.df$sample[samp.i])] = read.table(files.df$bc.qc.bg[samp.i], colClasses=c(rep("NULL",3),"numeric"))[,1]
        }
        write.table(bc.df, file=outfile.prefix, quote=FALSE, row.names=FALSE, sep="\t")
    } else {
        nb.chunks = ceiling(nrow(bin.df)/chunk.size)
        bin.df$chunk = rep(1:nb.chunks,each=chunk.size)[1:nrow(bin.df)]
        analyze.chunk <- function(df){
            ch.nb = as.numeric(df$chunk[1])
            bc.df = createEmptyDF(c("character",rep("integer",2), rep("numeric",nrow(files.df))), nrow(df))
            colnames(bc.df) = c("chr","start","end", as.character(files.df$sample))
            bc.df$chr = df$chr
            bc.df$start = df$start
            bc.df$end = df$end
            for(samp.i in 1:nrow(files.df)){
                bc.df[,as.character(files.df$sample[samp.i])] = read.bedix(files.df$bc.qc.bg[samp.i], df)[,4]
            }
            write.table(bc.df, file=outfile.prefix, quote=FALSE, row.names=FALSE, sep="\t", append=ch.nb>1, col.names=ch.nb==1)
            return(bc.df[sample(1:nrow(bc.df),chunk.size/nb.chunks),])
        }
        bc.df = dplyr::do(dplyr::group_by(bin.df,chunk),analyze.chunk(.))
    }
    
    if(appendIndex.outfile){
        Rsamtools::bgzip(outfile.prefix, overwrite=TRUE)
        file.remove(outfile.prefix)
        Rsamtools::indexTabix(paste(outfile.prefix,".bgz",sep=""), format="bed")
    }

    cor.bs <- function(x,nbsim=10,prop=.1,replace=FALSE){
        ni.sim = nrow(x)*prop
        cor.sim = sapply(1:nbsim,function(i)as.vector(cor(x[sample(1:nrow(x),ni.sim,replace=replace),],use="comp")))
        res = matrix(apply(cor.sim,1,median),ncol(x),ncol(x))
        colnames(res) = rownames(res) = colnames(x)
        return(res)
    }

    all.samples = colnames(bc.df)[-(1:3)]
    corbs = cor.bs(bc.df[,all.samples])
    cor.pw = corbs
    diag(cor.pw) = NA
    if(is.null(ref.samples)) ref.samples = colnames(cor.pw)
    meanCor = apply(cor.pw[ref.samples,ref.samples],1,function(e)mean(e,na.rm=TRUE))
    ## PCA
    pc = prcomp(t(na.exclude(bc.df[,all.samples])))
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
                pc.1.3=pc$x[,1:3]), cor.pw=corbs)
}
