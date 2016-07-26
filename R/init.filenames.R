##' Init the names for the files to use in the analysis to come. Specifically the
##' name of the files with the bin counts and compressed bin counts are created. Normalized bin counts, Z-scores, fold-changes will also be written following the paths created by this function.
##' @title Init file names for analysis
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample' and 'bam' are required.
##' @param dest.folder the folder where the files with the bin counts are to
##' be created.
##' @param sample.folder if TRUE (default), a folder is created for each sample
##' on the specified destination folder.
##' @param code if not NULL, a name for the analysis. Useful when PopSV is to be run
##' several times with different parameters. Default is NULL.
##' @param dest.folder.relative.path Is the specified folder defined relatively to
##' the working directory. Default is TRUE. Set to FALSE if providing absolute
##' paths (recommended).
##' @param sub.folder if non-null, sample folders will be structure into sub folder
##' defined by the column in 'files.df' named according to 'sub.folder's value. Useful when many samples
##' from different projects are studied together.
##' @return an updated data.frame with the information about the files.
##' For example a new column 'bc', 'bc.gz', 'bc.gc' and 'bc.gc.bz' holds the name
##' of the bin counts (raw/compressed), GC corrected bin counts (raw/compressed) files
##' to be created later. 
##' @author Jean Monlong
##' @export
init.filenames <- function(files.df, dest.folder = ".", sample.folder = TRUE, code = NULL, 
    dest.folder.relative.path = TRUE, sub.folder = NULL) {
    if (!all(c("sample", "bam") %in% colnames(files.df))) {
        stop("Columns 'sample' and 'bam' are required in 'files.df' data.frame.")
    }
    
    ## Convert potential factors into character Convert sample names into R-friendly
    ## names
    files.df$sample = make.names(as.character(files.df$sample))
    files.df$bam = as.character(files.df$bam)
    
    ## Check duplicate sample names
    files.df = unique(files.df)
    if (any(duplicated(files.df$sample))) {
        stop("Duplicated sample names: ", utils::head(files.df$sample[duplicated(files.df$sample)]))
    }
    
    ## Create folder structure as required
    if (dest.folder.relative.path) {
        dest.folder = paste(getwd(), dest.folder, sep = .Platform$file.sep)
    }
    sample.destf = dest.folder
    if (!is.null(sub.folder)) {
        if (all(!(sub.folder %in% colnames(files.df)))) {
            stop(paste("If non-null, 'sub.folder' must be the name of a column from 'files.df'. There is no column named", 
                sub.folder))
        }
        sample.destf = paste(sample.destf, files.df[, sub.folder], sep = .Platform$file.sep)
    }
    if (sample.folder) {
        sample.destf = paste(sample.destf, files.df$sample, sep = .Platform$file.sep)
    }
    ## Check if folder exists
    tmp = sapply(sample.destf, function(fold) {
        if (!file.exists(fold)) {
            dir.create(fold, recursive = TRUE)
        }
    })
    
    ## Create file names
    if (!is.null(code)) {
        code = paste0("-", code)
    }
    file.prefix = paste0(sample.destf, .Platform$file.sep, files.df$sample, code)
    files.df$bc = paste0(file.prefix, "-bc.tsv")
    files.df$bc.gz = paste0(files.df$bc, ".bgz")
    files.df$bc.gc = paste0(file.prefix, "-bc-gcCor.tsv")
    files.df$bc.gc.gz = paste0(files.df$bc.gc, ".bgz")
    files.df$bc.gc.norm = paste0(file.prefix, "-bc-gcCor-norm.tsv")
    files.df$bc.gc.norm.gz = paste0(files.df$bc.gc.norm, ".bgz")
    files.df$z = paste0(file.prefix, "-z.tsv")
    files.df$z.gz = paste0(files.df$z, ".bgz")
    files.df$fc = paste0(file.prefix, "-fc.tsv")
    files.df$fc.gz = paste0(files.df$fc, ".bgz")
    
    return(files.df)
} 
