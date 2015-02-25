##' Create Rscript file to run manually the bin counting for each sample.
##' If many samples and/or high coverage, and \code{BatchJobs} package is
##' not set up, this function quickly creates R scripts than can then be
##' send manually to a cluster.
##' @title Create Rscript files for manual bin counting
##' @param files.df a data.frame with the information about the files to
##' use. Columns 'sample', 'bam' and 'bc' are required and should be
##' present  after running 'initFileNames' function.
##' @param bin.df a data.frame with the information about the bins. Columns
##' 'chr', 'start' and 'end' are required.
##' @param dest.folder the folder where the R scripts files and other temporary
##' files are created.
##' @param proper if TRUE (default), reads with properly mapped pairs are
##' counted. If FALSE, reads with improper mapping are counted.
##' @param map.quality the minimum mapping quality (PHRED) for the reads
##' to count. Default is 30.
##' @return an updated data.frame with the information about the files.
##' Specifically a new column 'binBamRscripts' holds the name of the script
##' to run for each sample.
##' @author Jean Monlong
##' @export
createBinBamRscripts <- function(files.df, bin.df, dest.folder = ".", proper = TRUE, 
    map.quality = 30) {
    binf = paste(dest.folder, "binDef.RData", sep = .Platform$file.sep)
    save(bin.df, file = binf)
    create.script.samp <- function(samp, bamf, bcf) {
        scriptf = paste(dest.folder, .Platform$file.sep, "rscript-binBam-", samp, 
            ".R", sep = "")
        cat(paste0("library(PopSV)\nload(\"", binf, "\")\nbinBam(\"", bamf, "\", bin.df, \"", 
            bcf, "\", proper=", proper, ", map.quality=", map.quality, ")\n"), file = scriptf)
        return(scriptf)
    }
    files.df$binBamRscripts = sapply(1:nrow(files.df), function(samp.i) {
        create.script.samp(files.df$sample[samp.i], files.df$bam[samp.i], files.df$bc[samp.i])
    })
    return(files.df)
} 
