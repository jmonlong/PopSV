##' Computes the diversity of the supporting bins used during normalization. For now it only computes the number of bins in the same chromosome as the bin to normalize. To avoid over-normalization of large CNVs or full-chromosome aneuploidy (e.g. in cancer), we would like the majority of the supporting to be located different chromosomes than the bin to normalize.
##' @title Targeted normalization QC: diversity of the supporting bins
##' @param norm.stats the name of the file with the normalization statistics ('norm.stats' in 'tn.norm' function) or directly a 'norm.stats' data.frame.
##' @param out.pdf the name of the PDF file to create.
##' @param chunk.size the size of a chunk to read the file. Default is 1000.
##' @param nb.cores the number of cores to use. Default is 1.
##' @return a data.frame with the bins location and the computed metric.
##' @author Jean Monlong
##' @export
tn.norm.qc.div <- function(norm.stats, out.pdf = "normStats-QC-supportDiversity.pdf", chunk.size=1e3, nb.cores=1){

  read.chunk <- function(chunk.start=NULL, chunk.end=NULL){
    col.n = utils::read.table(norm.stats, nrows=1, sep="\t", header=FALSE, as.is=TRUE)
    dt = suppressWarnings(data.table::fread(norm.stats,nrows=chunk.end-chunk.start+1, skip=chunk.start, header=FALSE, sep="\t", showProgress=FALSE))
    data.table::setnames(dt, as.character(col.n))
    as.data.frame(dt)
  }

  ## Row number
  ns.1 = data.table::fread(norm.stats, header = TRUE, select=1, showProgress=FALSE)
  nrows = nrow(ns.1)
  rm(ns.1)

  ## Compute chunk index
  if(!is.null(chunk.size) && chunk.size<nrows){
    chunks = tapply(1:nrows, rep(1:ceiling(nrows/chunk.size), each=chunk.size)[1:nrows], identity)
  } else {
    chunks = list(1:nrows)
  }

  ## Analyze each chunk
  div.df = parallel::mclapply(1:length(chunks), function(ch.ii){
    ## Read chunk
    res.df = read.chunk(min(chunks[[ch.ii]]),max(chunks[[ch.ii]]))
    sup.ii = setdiff(colnames(res.df),  c("chr", "start", "end", "d.max", "m", "sd", "nb.remove"))
    binDiversity <- function(x){
      chr = x[1]
      chr.sup = unlist(lapply(strsplit(x[-1], "-"), "[", 1))
      return(sum(chr.sup==chr, na.rm=TRUE))
    }
    res.df$inter.chr = apply(res.df[, c("chr", sup.ii)], 1, binDiversity)
    res.df$nb.chr = apply(res.df[, sup.ii], 1, function(x)length(unique(x)))
    res.df[,c("chr","start","end","inter.chr","nb.chr","d.max", "m")]
  }, mc.cores=nb.cores)
  div.df = as.data.frame(data.table::rbindlist(div.df))

  ## Graph
  if(!is.null(out.pdf)){
      grDevices::pdf(out.pdf, 8, 6)
    chr = inter.chr = NULL ## Uglily appeases R checks
      print(ggplot2::ggplot(div.df[which(!is.na(div.df$d.max) & div.df$d.max>=0),], ggplot2::aes(x=chr, y=inter.chr)) + ggplot2::geom_violin(scale="width") + ggplot2::theme_bw() + ggplot2::ylab("number of supporting bins in same chr as normalized bin"))
        grDevices::dev.off()

  }

  return(div.df)
}
