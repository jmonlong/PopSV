##' Init the names for the files to use in the analysis to come. Specifically the
##' name of the files with the bin counts and compressed bin counts are created.
##' @title Init file names for analysis
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bam' are required.
##' @param dest.folder the folder where the files with the bin counts are to
##' be created.
##' @param sample.folder if TRUE (default), a folder is created for each sample
##' on the specified destination folder.
##' @param code if not NULL, a name for the analysis. Useful PopSV is to be run
##' several times with different parameters. Default is NULL.
##' @param dest.folder.relative.path Is the specified folder defined relatively to
##' the working directory. Default is TRUE. Set to FALSE if providing absolute
##' paths (recommended).
##' @return an updated data.frame with the information about the files.
##' Specifically a new column 'bc' and 'bc.gz' holds the name of the bin counts and
##' compressed bin counts files to be created later.
##' @author Jean Monlong
##' @export
initFileNames <- function(files.df, dest.folder=".", sample.folder=TRUE, code=NULL, dest.folder.relative.path=TRUE){
    if(!all(c("sample","bam")%in%colnames(files.df))){
        stop("Columns 'sample' and 'bam' are required in 'files.df' data.frame.")
    }
    files.df$sample = as.character(files.df$sample)
    files.df$bam = as.character(files.df$bam)
    if(dest.folder.relative.path){
        dest.folder = paste(getwd(),dest.folder,sep=.Platform$file.sep)
    }
    sample.destf = dest.folder
    if(sample.folder){
        sample.destf = paste(sample.destf,files.df$sample,sep=.Platform$file.sep)
    }
    ## Check if folder exists
    tmp = sapply(sample.destf, function(fold){
        if(!file.exists(fold)){
            dir.create(fold,recursive=TRUE)
        }
    })
    if(!is.null(CODE)){
        CODE = paste0("-",CODE)
    }
    files.df$bc = paste0(dest.folder,.Platform$file.sep,files.df$sample,CODE,"-bc.tsv")
    files.df$bc.gz = paste0(files.df$bc,".gz")
    return(files.df)
}
