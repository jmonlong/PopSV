##' Split join results into one file per sample per feature. For example to split and write the Z-scores of reference samples after joined computation.
##' @title Split and write results in one file per sample
##' @param res a list with the results to split/write.
##' @param files.df a data.frame with information about each output file.
##' @param samples the names of the samples to use.
##' @param files.col the name of 'files.df' column to use for each 'res' element.
##' @param compress.index should the output be compressed and indexed. Default is TRUE.
##' @param nb.cores the number of computing cores to use. Default is 1.
##' @param append should the results be appended at the end of existing files. Default is FALSE.
##' @param reorder If 'compress.index=TRUE', should the files be reordered before compression/indexing. Default is FALSE.
##' @return 'Done' if everything worked fine.
##' @author Jean Monlong
##' @export
write.split.samples <- function(res, files.df, samples=NULL, files.col = c("z","fc"), compress.index = TRUE, nb.cores=1, append=FALSE, reorder=FALSE) {
  if (length(res) != length(files.col)) {
    stop("'res' and 'files.col' have different length.")
  }
  if(length(res[[1]])==1 & is.character(res[[1]])){
    
    read.chunk <- function(fn, chunk.start=NULL, chunk.end=NULL, samples=NULL){
      col.n = utils::read.table(fn, nrows=1, sep="\t", header=FALSE, as.is=TRUE)
      if(is.null(samples)){
        col.ii = 1:length(col.n)
      } else {
        col.ii = which(as.character(col.n) %in% c("chr","start","end",samples))
      }
      col.n = col.n[,col.ii]
      dt = suppressWarnings(data.table::fread(fn,nrows=chunk.end-chunk.start+1, skip=chunk.start, header=FALSE, sep="\t", select=col.ii))
      data.table::setnames(dt, as.character(col.n))
      as.data.frame(dt)
    }
    
    ## row number
    bc.1 = data.table::fread(res[[1]], header = TRUE, select=1)
    nrows = nrow(bc.1)
    rm(bc.1)

    ## Compute chunk index
    chunk.size = 1e4
    if(!is.null(chunk.size) && chunk.size<nrows){
      chunks = tapply(1:nrows, rep(1:ceiling(nrows/chunk.size), each=chunk.size)[1:nrows], identity)
    } else {
      chunks = list(1:nrows)
    }

    tmp = lapply(1:length(chunks), function(ch.ii){
      chunk.l = lapply(1:length(files.col), function(f.ii){
        read.chunk(res[[f.ii]], min(chunks[[ch.ii]]),max(chunks[[ch.ii]]), samples)
      })
      write.split.samples(chunk.l, files.df, samples=samples, files.col=files.col, compress.index=FALSE, append=append | ch.ii>1)
    })
    
  } else {
    if(is.null(samples)){
      samples = setdiff(colnames(res[[1]]), c("chr","start","end"))
    }
    tmp = parallel::mclapply(samples, function(samp) {
      for (ii in 1:length(files.col)) {
        res.f = res[[ii]][, c("chr", "start", "end", samp)]
        colnames(res.f)[4] = files.col[ii]
        res.f = with(res.f, dplyr::arrange(res.f, chr, start))
        utils::write.table(res.f, file = files.df[which(files.df$sample == samp), files.col[ii]], row.names = FALSE, quote = FALSE, sep = "\t", append=append, col.names=!append)
      }
    }, mc.cores=nb.cores)
  }

  if (compress.index) {
    files.tc = as.character(unlist(files.df[which(files.df$sample %in% samples), files.col]))
    comp.index.files(files.tc, reorder=reorder)
  }
  
  return("Done")
}
