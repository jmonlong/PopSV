##' Split join results into one file per sample per feature. For example to split and write the Z-scores of reference samples after joined computation.
##' @title Split and write results in one file per sample
##' @param res a list with the results to split/write.
##' @param files.df a data.frame with information about each output file.
##' @param samples the names of the samples to use.
##' @param files.col the name of 'files.df' column to use for each 'res' element.
##' @param compress.index should the output be compressed and indexed. Default is TRUE.
##' @param nb.cores the number of computing cores to use. Default is 1.
##' @param append should the results be appended at the end of existing files. Default is FALSE.
##' @return 'Done' if everything worked fine.
##' @author Jean Monlong
##' @export
write.split.samples <- function(res, files.df, samples, files.col = c("z","fc"), compress.index = TRUE, nb.cores=1, append=FALSE) {
  if (length(res) != length(files.col)) {
    stop("'res' and 'files.col' have different length.")
  }

  tmp = parallel::mclapply(samples, function(samp) {
    for (ii in 1:length(files.col)) {
      res.f = res[[ii]][, c("chr", "start", "end", samp)]
      colnames(res.f)[4] = files.col[ii]
      res.f = with(res.f, dplyr::arrange(res.f, chr, start))
      write.table(res.f, file = files.df[which(files.df$sample == samp), files.col[ii]], row.names = FALSE, quote = FALSE, sep = "\t", append=append, col.names=!append)
    }
  }, mc.cores=nb.cores)

  if (compress.index) {
    files.tc = as.character(unlist(files.df[which(files.df$sample %in% samples),
      files.col]))
    comp.index.files(files.tc)
  }

  return("Done")
}
