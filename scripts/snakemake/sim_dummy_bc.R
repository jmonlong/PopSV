library(dplyr)

## Bin counts for samples 2:20
S = 1000000
W = 500
sampids = 2:20
temp = lapply(sampids, function(samp){
  bc.df = data.frame(chr=rep(1:2, each=S/W), start=rep(seq(1, S, W), 2))
  bc.df$end = bc.df$start + W - 1
  bc.df$bc = rnorm(nrow(bc.df), 500, 100)
  bc.df$bc = ifelse(bc.df$bc<0, 0, bc.df$bc)
  out.con = gzfile(paste0('bc-samp', samp, '.tsv.gz'), 'wb')
  write.table(bc.df, file=out.con, quote=FALSE,
              sep='\t', row.names=FALSE)
  close(out.con)
})

## Sample info
info = tibble(sample=paste0('samp', 1:20), bam=NA, bc=paste0('bc-samp', 1:20, '.tsv.gz'))
info$bam[1] = 'samp1.bam'
info$bc[1] = NA
info$reference = TRUE
info$reference[19:20] = FALSE
write.table(info, file='info.tsv', sep='\t', quote=FALSE, row.names=FALSE)
