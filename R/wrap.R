##' Meant to be used from the command line with 'Rscript "PopSV::wrap()" countsample samp1 config.RData bins.RData'.
##' @title Wrapper used for making pipelines
##' @param args command line arguments
##' @return the status of the computation.
##' @author Jean Monlong
##' @export
wrap <- function(args=commandArgs(TRUE)){

  ## Create the config file from a TSV file with sample name and bam paths
  initfilename <- function(args){
    sample_bam_file = args[1]
    popsv_config_file = args[2]
    ref_samps_file = ifelse(length(args)>2, args[3], NA)
    message('Create configuration file and folders')
    bam.files = utils::read.table(sample_bam_file, as.is=TRUE, header=TRUE)
    files.df = init.filenames(bam.files)
    if(!is.na(ref_samps_file) && file.exists(ref_samps_file)){
      message('Read the list of reference sample names')
      ref.samps = make.names(scan(ref_samps_file, ''))
      files.df$reference = files.df$sample %in% ref.samps
    }    
    save(files.df, file=popsv_config_file)
    return('Done')
  }

  ## Prepare the bins.
  ## Add a column with GC content to the bins.
  ## If the bin file doesn't exist, bins with 'bin_size' will be generated.
  ## 'genome' either 'hg19' or 'GRCh38'
  prepbins <- function(args){
    bin_file_input = args[1]
    bin_size = as.integer(args[2])
    nb_chunks = as.integer(args[3])
    bin_file = args[4]
    genome = ifelse(length(args)>4, args[5], 'hg19')
    message('Prepare bin file')
    if(file.exists(bin_file_input)){
      message('File exists. Reading ', bin_file_input)
      bins.df = utils::read.table(bin_file_input, as.is=TRUE, header=TRUE, sep='\t')
      bins.df$start = as.integer(bins.df$start)
      bins.df$end = as.integer(bins.df$end)
    } else {
        message('Binning in windows of ', bin_size, ' bp')
        bins.df = fragment.genome(bin_size, genome=genome)
    }
    if(!('GCcontent' %in% colnames(bins.df))){
      message('Retrieve GC content')
      bins.df = getGC(bins.df, genome=genome)
    }
    message('Prepare chunks of bins')
    sm.chunk.size = ceiling(nrow(bins.df)/nb_chunks)
    bins.df = chunk.bin(bins.df, bg.chunk.size=5e5, sm.chunk.size=sm.chunk.size)
    save(bins.df, file=bin_file)
    return('Done')
  }

  ## Count reads from a BAM file for a sample
  countsample <- function(args){
    bam.f = args[1]
    bins.df = NULL
    load(args[2])
    bc.f = args[3]
    bc.f = gsub('\\.bgz$', '', bc.f)
    bb.o = bin.bam(bam.f, bins.df, bc.f)
  }

  ## Correct the read counts for GC bias for a sample
  gccorrect <- function(args){
    bins.df = NULL
    load(args[2])
    outfile.nogz = gsub('\\.bgz$', '', args[3])
    correct.GC(args[1], bins.df, outfile.nogz)
    return('Done')
  }

  ## Merge read counts for reference samples
  preprefs <- function(args){
    popsv_config_file = args[1]
    bin_file = args[2]
    ref_file = args[3]
    cont_sample_file = args[4]
    nb.cores = ifelse(length(args)>4, args[5], 1)
    graph_out = 'preprefs.pdf'
    max_nb_refs=200
    bins.df = files.df = NULL
    ref_file = gsub('\\.bgz$', '', ref_file)
    bins.df = NULL
    load(bin_file)
    files.df = utils::read.table(popsv_config_file, as.is=TRUE, header=TRUE, sep='\t')
    files.df = files.df[which(as.logical(files.df$reference)),]
    grDevices::pdf(graph_out)
    qc.o = qc.samples(files.df, bins.df, ref_file, nb.ref.samples=max_nb_refs,
                      nb.cores=nb.cores)
    grDevices::dev.off()
    write(qc.o$cont.sample, file=cont_sample_file)
  }

  ## Normalize the reference samples
  normrefs <- function(args){
    ref_file = args[1]
    bin_file = args[2]
    cont_sample_file = args[3]
    chunk = as.integer(args[4])
    res_file = args[5]
    nb.support.bins = as.integer(args[6])
    bins.df = NULL
    load(bin_file)
    chunks = sort(unique(as.character(bins.df$sm.chunk)))
    chunk = chunks[chunk]
    bins.df.chunk = bins.df[which(bins.df$sm.chunk==chunk),]
    bg.chunk = utils::head(bins.df.chunk$bg.chunk, 1)
    bc.df = read.bedix(ref_file, bins.df[which(bins.df$bg.chunk==bg.chunk),])
    cont.sample = scan(cont_sample_file, '')
    res = tn.norm(bc.df, cont.sample, bins=bins.df.chunk$bin,
                  nb.support.bins=nb.support.bins)
    utils::write.table(res$norm.stats, file=res_file, sep='\t', row.names=FALSE, quote=FALSE)
  }

  ## ## Merge the output of the normalization step
  ## mergeoutrefs <- function(args){
  ##   norm_ref_prefix = args[1]
  ##   nb_chunks = as.integer(args[2])
  ##   ref_prefix = args[3]
  ##   done = ifelse(length(args)>3, args[4], NA)
  ##   norm_ref_files = paste0(norm_ref_prefix, '_', 1:nb_chunks, '.RData')
  ##   outfile = paste0(ref_prefix, 'norm-stats.tsv')
  ##   if(file.exists(outfile)){
  ##     file.remove(outfile)
  ##   }
  ##   tmp = lapply(norm_ref_files, function(ff){
  ##     res = NULL
  ##     load(ff)
  ##     utils::write.table(res$norm.stats, file=outfile, sep='\t', row.names=FALSE,
  ##                        append=file.exists(outfile), col.names=!file.exists(outfile),
  ##                        quote=FALSE)
  ##   })
  ##   file.remove(norm_ref_files)
  ##   if(!is.na(done)) cat("Done", file=done)
  ##   message('Done')
  ## }

  ## Call CNV in a sample
  callsample <- function(args){
    samp_name = args[1]
    popsv_config_file = args[2]
    bin_file = args[3]
    cont_sample_file = args[4]
    ref_file = args[5]
    norm_stats = args[6]
    FDR_th = as.numeric(args[7])
    cnv_file = args[8]
    samp_name = make.names(samp_name)
    bins.df = files.df = NULL
    files.df = utils::read.table(popsv_config_file, as.is=TRUE, header=TRUE, sep='\t')
    bins.df = NULL
    load(bin_file)
    cont.sample = scan(cont_sample_file, '')
    if(!file.exists(files.df$z[which(files.df$sample == samp_name)])){
      if(!file.exists(ref_file)){
        message('Adding bgz')
        ref_file = paste0(ref_file, '.bgz')
      }
      message('Normalize bin count and compute Z-score')
      tn.test.sample(samp_name, files.df, cont.sample, ref_file,
                     norm_stats, z.poisson=TRUE,
                     aberrant.cases=FALSE)
    }
    message('Call CNVs')
    stitch.dist = mean(bins.df$end-bins.df$start+1)*2
    sub.z = min(1e3, nrow(bins.df)/3)
    cnv.df = call.abnormal.cov(files.df=files.df, samp=samp_name, FDR.th=FDR_th,
                               merge.cons.bins="cbs", z.th="sdest",
                               norm.stats=norm_stats,
                               stitch.dist=stitch.dist, gc.df=bins.df,
                               min.normal.prop=.6, sub.z=sub.z)
    message('Write CNV output')
    utils::write.table(cnv.df, file=cnv_file, quote=FALSE, sep='\t', row.names=FALSE)
    message('Done')
  }

  ## Read arguments and call function
  if(args[1] == 'initfilename'){
    initfilename(args[-1])
  } else if(args[1] == 'prepbins'){
    prepbins(args[-1])
  } else if(args[1] == 'countsample'){
    countsample(args[-1])
  } else if(args[1] == 'gccorrect'){
    gccorrect(args[-1])
  } else if(args[1] == 'preprefs'){
    preprefs(args[-1])
  } else if(args[1] == 'normrefs'){
    normrefs(args[-1])
  ## } else if(args[1] == 'mergeoutrefs'){
  ##   mergeoutrefs(args[-1])
  } else if(args[1] == 'callsample'){
    callsample(args[-1])
  } else {
    stop('Unknown step: ', args[1])
  }

  return('Done')
}
