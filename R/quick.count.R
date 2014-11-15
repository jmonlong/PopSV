##' This function retrieves bin count for all samples in the input file information
##' data.frame. It was designed to count reads in a small amount of bins (a subset
##' of what will be used for the full analysis). This counts can be used later to
##' compare a potentially large number samples and guide the selection of samples
##' to consider in the analysis. 
##' @title Counts reads across samples in a small number of bins
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bam' are required.
##' @param bins.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @return a data.frame with the bin location and read counts for all samples.
##' @author Jean Monlong
##' @export
quick.count <- function(files.df, bins.df, nb.cores=1){
    if(nb.cores>1){
        bc.l = parallel::mclapply(files.df$sample, function(samp){
            bc.s = bin.bam(files.df$bam[files.df$sample==samp], bins.df, appendIndex.outfile=FALSE, chunk.size=nrow(bins.df))
            bc.s$bc
        },mc.cores=nb.cores)
    } else {
        bc.l = lapply(files.df$sample, function(samp){
            bc.s = bin.bam(files.df$bam[files.df$sample==samp], bins.df, appendIndex.outfile=FALSE, chunk.size=nrow(bins.df))
            bc.s$bc
        })
    }
    bc.df = createEmptyDF(c(sapply(c("chr","start","end"), function(cn)class(bins.df[,cn])), rep("integer",nrow(files.df))),nrow(bins.df))
    colnames(bc.df) = c("chr","start","end",as.character(files.df$sample))
    bc.df[,1:3] = bins.df[,c("chr","start","end")]
    for(ii in 1:length(bc.l)){
        bc.df[,3+ii] = bc.l[[ii]]
    }
    bc.df
}
