##' Normalize bin counts.
##' @title Normalize bin counts
##' @param bc a data.frame with the bin counts. 
##' @param files.df a data.frame with the names of the files with the the bin counts. Used if 'bc' is NULL.
##' @param bins.df a data.frame with the information about the bins. Columns 'chr', 'start'
##' and 'end' are required.
##' @param method the normalization method: 'tn' for targeted normalization; 'pca' for
##' principal component regression; 'medvar' for median/variance global normalization.
##' @param ... extra parameters for the normalization methods.
##' @return the output from 'tn.norm'/'pca.norm'/'medvar' function, depending of
##' the choice of 'method'.
##' @author Jean Monlong
##' @export
norm.bc <- function(bc=NULL, files.df=NULL, bins.df=NULL, method=c("tn","pca","medvar"), ...){
    if(is.null(bc)){
        if(is.null(files.df)){
            stop("Please provide 'bc' or 'files.df' (with eventually 'bins.df')")
        }

        head.df = read.table(files.df$bc.gc.gz[1], as.is=TRUE, nrows=10, header=TRUE)
        types = apply(head.df, 2, typeof)

        if(is.null(bins.df)){
            bc.s = read.table(files.df$bc.gc.gz[1], as.is=TRUE, header=TRUE)
            samp.int = 2:nrow(files.df)
            nb.bins = nrow(bc.s)
        } else {
            samp.int = 1:nrow(files.df)
            nb.bins = nrow(bins.df)
        }

        bc = createEmptyDF(c(types[-length(types)],rep("numeric", nrow(files.df))), nb.bins)
        colnames(bc) = c(colnames(head.df)[-length(types)], files.df$sample)
        
        if(is.null(bins.df)){
            bc[,as.character(files.df$sample[1])] = bc.s$bc
        }
        
        for(samp.i in samp.int){
            if(is.null(bins.df)){
                bc.s = read.table(files.df$bc.gc.gz[samp.i], header=TRUE, as.is=TRUE)
            } else {
                bc.s = read.bedix(files.df$bc.gc.gz[samp.i], bins.df)
            }
            bc[,as.character(files.df$sample[samp.i])] = bc.s$bc
        }
    }

    if(method[1]=="tn"){
        return(tn.norm(bc,...))
    } else if(method[1]=="pca"){
        return(pca.norm(bc,...))
    } else if(method[1]=="medvar"){
        return(medvar.norm(bc,...))
    } else {
        stop("Please choose a method from 'tn', 'pca' or 'medvar'")
    }
}

