##' Split join results into one file per sample per feature. For example to split and write the Z-scores of reference samples after joined computation.
##' @title Split and write results in one file per sample
##' @param res a list with the results to split/write.
##' @param files.df a data.frame with information about each output file.
##' @param samples the names of the samples to use.
##' @param res.n the name of 'res' elements containing joined results.
##' @param files.col the name of 'files.df' to use for each 'res' element.
##' @param compress.index should the output be compressed and indexed. Default is TRUE.
##' @return "Done" if everything worked fine. 
##' @author Jean Monlong
##' @export
write.split.samples <- function(res, files.df, samples, res.n=c("z","fc"), files.col=c("z","fc"),compress.index=TRUE){

    if(length(res.n)!=length(files.col)) stop("'res.n' and 'files.col' have different length.")
    
    tmp = sapply(samples, function(samp){
        for(ii in 1:length(files.col)){
            res.f = res[[res.n[ii]]][,c("chr","start","end",samp)]
            colnames(res.f)[4] = res.n[ii]
            res.f = with(res.f, dplyr::arrange(res.f, chr, start))
            write.table(res.f, file=files.df[which(files.df$sample==samp),files.col[ii]], row.names=FALSE, quote=FALSE,sep="\t")
        }
    })

    if(compress.index){
        files.tc = as.character(unlist(files.df[which(files.df$sample %in% samples),files.col]))
        comp.index.files(files.tc)
    }
    
    return("Done")
}
