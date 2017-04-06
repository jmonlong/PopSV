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
    chunk.size = 1e4
    if(length(res[[1]])==1 & is.character(res[[1]])){

        con.l = lapply(res, file, "r")
        header.l = lapply(con.l, function(con){
                              unlist(strsplit(readLines(con, n = 1), "\t"))
                          })
        if(is.null(samples)){
            samples = setdiff(header.l[[1]], c("chr","start","end"))
        }

        read.chunk <- function(){
            res.l = lapply(1:length(con.l), function(con.ii){
                               res = tryCatch(utils::read.table(con.l[[con.ii]], as.is=TRUE, nrows=chunk.size), error=function(e)return(NULL))
                               if(is.null(res)){
                                   return(NULL)
                               }
                               colnames(res) = header.l[[con.ii]]
                               col.ii = which(colnames(res) %in% c("chr","start","end",samples))
                               res[,col.ii]
                           })
            if(is.null(res.l[[1]])){
                return(NULL)
            }
            return(res.l)
        }

        firstChunk = TRUE
        while(!is.null((chunk.o = read.chunk()))){
            write.split.samples(chunk.o, files.df, samples=samples, files.col=files.col, compress.index=FALSE, append=append | !firstChunk)
            firstChunk = FALSE
        }

        close.con = lapply(con.l, close)
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
