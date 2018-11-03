##' Transforms the BED/data.frame format into a VCF format and write a VCF file.
##' The calls are merged if they have exactly the same coordinates ('chr', 'start' and 'end').
##' @title Write VCF file
##' @param cnv.df the data.frame with the cnv calls (from 'call.abnormal.cov').
##' @param output.file the file name for the output VCF file. Default is 'cnvs.vcf'
##' @param genome BSgenome to retrieve REF nucleotide. Default is hg19 (BSgenome.Hsapiens.UCSC.hg19::Hsapiens).
##' @return the output file name.
##' @author Jean Monlong
writeVcf <- function(cnv.df, output.file='cnvs.vcf', genome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){
  ## Checking input
  if(!all(c('chr', 'start', 'end', 'nb.bin.cons', 'fc', 'sample') %in% colnames(cnv.df))){
    stop("Missing column. 'cnv.df' must have 'chr', 'start', 'end', 'nb.bin.cons', 'sample' and 'fc' columns")
  }

  ## Meta infos and headers
  meta = c("##fileformat=VCFv4.2",
    paste0("##fileDate=", format(Sys.time(), "%Y%m%d")),
    '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">',
    '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##ALT=<ID=CNV,Description="Copy number variable region">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype">',
    '##FORMAT=<ID=CNE,Number=1,Type=Float,Description="Copy number estimate">')
  headers = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

  ## Estimate bin size for the breakpoint confidence intervals
  b1.ii = which(cnv.df$nb.bin.cons==1)
  binsize = min(cnv.df$end[b1.ii] - cnv.df$start[b1.ii])

  ## To find nucleotide for each position (for REF column)
  getRefBase <- function(df){
    chrs = df$chr
    if (!grepl("chr", chrs[1])) {
      chrs = paste("chr", chrs, sep = "")
    }
    nuc = Biostrings::getSeq(genome, chrs, df$start, df$start)
    df$ref = as.character(nuc)
    df
  }

  ## Create CNV genotypes
  cnv.df$cn = round(cnv.df$fc*2)
  cnv.df$cn = ifelse(cnv.df$cn==2, 2+sign(cnv.df$fc-1), cnv.df$cn)
  cnv.df$gt = paste0("./.:", cnv.df$cn, ':', round(cnv.df$fc*2,2))

  ## Merge samples
  gt = NULL ## Uglily appeases R checks
  cnv.m = cnv.df[,c("chr", "start", "end", "sample", "gt")]
  cnv.m = tidyr::spread(cnv.m, sample, gt, fill='./.:2:2')

  ## Update headers with sample names
  samples = setdiff(colnames(cnv.m), c("chr", "start", "end"))
  headers = c(headers, samples)
  headers = paste(headers, collapse="\t")

  ## VCF format data.frae
  vcf.df = getRefBase(cnv.m)
  vcf.df$id = '.'
  vcf.df$alt = '<CNV>'
  vcf.df$qual = '.'
  vcf.df$filter = 'PASS'
  vcf.df$info = paste0("SVTYPE=CNV;END=",vcf.df$end,";SVLEN=",
      vcf.df$end-vcf.df$start+1,";CIPOS=-",binsize/2,",", binsize/2,";CIEND=-",
      binsize/2,",", binsize/2)
  vcf.df$format='GT:CN:CNE'
  vcf.df$end = NULL
  vcf.df = vcf.df[,c("chr", "start", "id", "ref", "alt", "qual", "filter", "info", "format", samples)]

  ## Write file
  write(meta, file=output.file)
  write(headers, file=output.file, append=TRUE)
  utils::write.table(vcf.df, file=output.file, append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

  return(output.file)
}
