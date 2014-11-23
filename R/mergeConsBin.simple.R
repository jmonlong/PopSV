##' Merge consecutive bins (end1+1->start2).
##' @title Simple merge of abnormal bins
##' @param res.df a data.frame with the Z-scores. Columns 'chr', 'start', 'end' and 'z'
##' are required.
##' @param col.mean a list of potential columns for which the mean value should be kept
##' when merging bins.
##' @return a data.frame similar to the input 'res.df' but with an extra 'nbBinCons'
##' column (the number of bin merged for each event).
##' @author Jean Monlong
##' @keywords internal
mergeConsBin.simple <- function(res.df,col.mean=c("z","pv","qv","fc")){

    link.annotate.f <- function(df){
        chr.f = df$chr[1]
        df = dplyr::arrange(df, start)
        cons.v = df$end[-nrow(df)] + 1 == df$start[-1]
        link.v = paste0(chr.f,"none",1:nrow(df))
        rl = rle(cons.v)
        rlv = rl$values
        cs.i = c(0,cumsum(rl$lengths))
        for(cons.i in which(rlv)){
            link.v[(cs.i[cons.i]+1):(cs.i[cons.i+1]+1)] = paste0(chr.f,"cons",cons.i)
        }
        df$link = as.character(link.v)
        return(df)
    }
    
    res.df = dplyr::do(dplyr::group_by(res.df,chr),link.annotate.f(.))
    merge.bin.f <- function(df){
        if(nrow(df)<3) return(data.frame())
        res = data.frame(chr=df$chr[1],
            start=df$start[1],
            end=df$end[nrow(df)],
            nbBinCons=nrow(df), stringsAsFactors=FALSE)
        cbind(res,t(apply(df[,intersect(colnames(df),col.mean),drop=FALSE],2,mean)))
    }
    res.df = dplyr::do(dplyr::group_by(res.df,link),merge.bin.f(.))
    res.df$link = NULL
    return(res.df)   
}
