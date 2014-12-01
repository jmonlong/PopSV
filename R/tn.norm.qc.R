##' Following targeted normalization, some Quality Control graphs are produced to see if some bins were problematic. TODO MORE
##' @title QC graphs for Targeted Normalization
##' @param norm.stats the name of the file with the normalization statistics ('norm.stats' in 'tn.norm' function) or directly a 'norm.stats' data.frame.
##' @param out.pdf the name of the PDF file to create.
##' @return the name of the created PDF file.
##' @author Jean Monlong
##' @export
tn.norm.qc <- function(norm.stats, out.pdf="normStats-QC.pdf"){
    ## load norm statistics
    if(is.character(norm.stats) & length(norm.stats)==1){
        headers = read.table(norm.stats,nrows=1,as.is=TRUE)
        colC = rep("NULL",length(headers))
        names(colC) = headers
        colC[c("chr","start","end","d.max","m","sd","nb.remove")] = c("character", rep("integer",2),rep("numeric", 4))
        res.df = read.table(norm.stats,header=TRUE,colClasses=colC)
    } else {
        res.df = norm.stats[,c("chr","start","end","d.max","m","sd","nb.remove")]
        rm(norm.stats)
    }

    res.df = subset(res.df, d.max!=-1 & !is.na(d.max))
    
    pdf(out.pdf, 8,6)
    print(ggplot2::ggplot(res.df, ggplot2::aes(x=d.max)) + ggplot2::geom_histogram() + ggplot2::theme_bw() + ggplot2::xlab("correlation distance to last supporting bin") + ggplot2::ylab("number of bins"))
    print(ggplot2::ggplot(res.df, ggplot2::aes(x=m)) + ggplot2::geom_histogram() + ggplot2::theme_bw() + ggplot2::xlab("average normalized coverage") + ggplot2::ylab("number of bins"))
    print(ggplot2::ggplot(res.df, ggplot2::aes(x=m+1)) + ggplot2::geom_histogram() + ggplot2::theme_bw() + ggplot2::scale_x_log10() + ggplot2::xlab("average normalized coverage") + ggplot2::ylab("number of bins"))
    print(ggplot2::ggplot(res.df, ggplot2::aes(x=nb.remove)) + ggplot2::geom_histogram() + ggplot2::theme_bw() + ggplot2::xlab("number of outlier samples removed") + ggplot2::ylab("number of bins"))
    ## Pairwise graphs
   print(ggplot2::ggplot(res.df, ggplot2::aes(x=d.max, y=m+1, fill=log10(..count..))) + ggplot2::stat_bin2d() + ggplot2::theme_bw() + ggplot2::xlab("correlation distance to last supporting bin") + ggplot2::ylab("average normalized coverage") + ggplot2::scale_y_log10() + ggplot2::scale_fill_gradient(name="log10(nb bins)",low="white", high="red"))
   print(ggplot2::ggplot(res.df, ggplot2::aes(x=d.max, y=nb.remove, fill=log10(..count..))) + ggplot2::stat_bin2d() + ggplot2::theme_bw() + ggplot2::xlab("correlation distance to last supporting bin") + ggplot2::ylab("number of outlier samples removed")  + ggplot2::scale_fill_gradient(name="log10(nb bins)",low="white", high="red"))
    print(ggplot2::ggplot(res.df, ggplot2::aes(x=m+1, y=nb.remove, fill=log10(..count..))) + ggplot2::stat_bin2d() + ggplot2::theme_bw() + ggplot2::xlab("average normalized coverage") + ggplot2::ylab("number of outlier samples removed")  + ggplot2::scale_fill_gradient(name="log10(nb bins)",low="white", high="red") + ggplot2::scale_x_log10())
    dev.off()
    
    return(out.pdf)
}
