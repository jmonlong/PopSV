message(" /!\ Still need to be tested /!\ ")

message("Two functions :
- 'autoGCcounts' to count BC in each sample.
- 'autoNormTest' to normalize and test all the samples.
")

autoGCcounts <- function(files.f, bins.f, redo=NULL, sleep=180, status=FALSE, file.suffix="", lib.loc=NULL, skip=NULL){
  load(files.f)
  
  message("\n== 1) Get bin counts in each sample and correct for GC bias.\n")
  stepName = paste0("getBC",file.suffix)
  if(any(redo==1)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(!any(skip==1) & length(findJobs(reg))==0){
    getBC.f <- function(file.i, bins.f, files.df){
      library(PopSV, lib.loc=lib.loc)
      load(bins.f)
      system(paste0("/icgc/icgc-storage-client/bin/icgc-storage-client --profile collab download --object-id ",files.df$object[file.i]," --output-dir /root/bams"))
      bam.f = paste0("/root/bams/",files.df$bam[file.i],"/",files.df$object[file.i])
      bai.f = paste0("/root/bams/",files.df$bam[file.i],".bai")
      bb.o = bin.bam(bam.f, bins.df, files.df$bc[file.i], bai.file=bai.f)
      correct.GC(files.df$bc.gz[file.i], bins.df, files.df$bc.gc[file.i])
      ## system(paste0("rm /root/bams/",files.df$bam[file.i]))
      bb.o
    }
    batchMap(reg, getBC.f,1:nrow(files.df), more.args=list(bins.f=bins.f, files.df=files.df))
    submitJobs(reg, findJobs(reg))
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  if(status) showStatus(reg)
}
  
autoNormTest <- function(files.f, bins.f, redo=NULL, rewrite=FALSE, sleep=180, status=FALSE, loose=FALSE, file.suffix="", lib.loc=NULL, ref.samples=NULL){
  load(files.f)
  
  message("\n== 1) Sample QC and reference definition.\n")
  bc.ref.f = paste0("bc-gcCor",file.suffix,".tsv")
  sampQC.pdf.f = paste0("sampQC",file.suffix,".pdf")
  stepName = paste0("sampQC",file.suffix)
  if(any(redo==1)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(length(findJobs(reg))==0){
    if(!is.null(ref.samples)){
      files.ref = subset(files.df, sample %in% ref.samples)
    } else {
      files.ref = files.df
    }
    sampQC.f <- function(bc.all.f, bins.f, files.df, sampQC.pdf.f, lib.loc){
      load(bins.f)
      library(PopSV, lib.loc=lib.loc)
      pdf(sampQC.pdf.f)
      qc.o = qc.samples(files.df, bins.df, bc.all.f, nb.cores=6, nb.ref.samples=200)
      dev.off()
      qc.o 
    }
    batchMap(reg, sampQC.f,bc.ref.f, more.args=list(bins.f=bins.f, files.df=files.ref, sampQC.pdf.f=sampQC.pdf.f, lib.loc=lib.loc))
    submitJobs(reg, 1)
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  samp.qc.o = loadResult(reg, 1)
  save(samp.qc.o, file=paste0(stepName,".RData"))
  if(status) showStatus(reg)

  message("\n== 2) Reference sample normalization.\n")
  stepName = paste0("bcNormTN",file.suffix)
  if(any(redo==2)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(length(findJobs(reg))==0){
    load(bins.f)
    if(all(colnames(bins.df)!="sm.chunk")){
      bins.df = chunk.bin(bins.df, bg.chunk.size=1e5, sm.chunk.size=1e4, large.chr.chunks=TRUE)
      save(bins.df, file=bins.f)
    }
    bcNormTN.f <- function(chunk.id, file.bc, file.bin, cont.sample, lib.loc){
      load(file.bin)
      library(PopSV, lib.loc=lib.loc)
      bc.df = read.table(file.bc, header=TRUE, as.is=TRUE)
      tn.norm(bc.df, cont.sample, bins=subset(bins.df, sm.chunk==chunk.id)$bin)
    }
    batchMap(reg, bcNormTN.f,unique(bins.df$sm.chunk), more.args=list(file.bc=samp.qc.o$bc, file.bin=bins.f,cont.sample=samp.qc.o$cont.sample, lib.loc=lib.loc))
    submitJobs(reg, findJobs(reg))
    waitForJobs(reg, sleep=sleep)
  } 
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  ## Write normalized bin counts and reference metrics
  out.files = paste(paste0("ref",file.suffix), c("bc-norm.tsv", "msd.tsv"), sep="-")
  if(rewrite | all(!file.exists(out.files))){
    if(any(file.exists(out.files))){
      tmp = file.remove(out.files[which(file.exists(out.files))])
    }
    tmp = reduceResultsList(reg, fun=function(res, job){
      write.table(res$bc.norm, file=out.files[1], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[1]), col.names=!file.exists(out.files[1]))
      write.table(res$norm.stats, file=out.files[2], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[2]), col.names=!file.exists(out.files[2]))
    })
  }
  if(status) showStatus(reg)
  
  message("\n== 3) Compute Z-scores in reference samples.\n")
  stepName = paste0("zRef",file.suffix)
  if(any(redo==3)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(length(findJobs(reg))==0){
    zRef.f <- function(bc.f, files.df, lib.loc){
      library(PopSV, lib.loc=lib.loc)
      z.comp(bc.f=bc.f, files.df=files.df, nb.cores=3, z.poisson=TRUE, chunk.size=1e3)
    }
    batchMap(reg, zRef.f,out.files[1], more.args=list(files.df=files.df, lib.loc=lib.loc))
    submitJobs(reg, 1)
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  if(status) showStatus(reg)

  message("\n== 4) Normalization and Z-score computation for other samples.\n")
  stepName = paste0("zOthers",file.suffix)
  if(any(redo==4)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(length(findJobs(reg))==0){
    callOthers.f <- function(samp, cont.sample, files.df, norm.stats.f, bc.ref.f, lib.loc){
      library(PopSV, lib.loc=lib.loc)
      tn.test.sample(samp, files.df, cont.sample, bc.ref.f, norm.stats.f, z.poisson=TRUE, aberrant.cases=TRUE)
    }
    batchMap(reg, callOthers.f,setdiff(files.df$sample, samp.qc.o$ref.samples), more.args=list(cont.sample=samp.qc.o$cont.sample, files.df=files.df, norm.stats.f=out.files[2], bc.ref.f=samp.qc.o$bc, lib.loc=lib.loc))
    submitJobs(reg, findJobs(reg))
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  if(status) showStatus(reg)

  message("\n== 5) Calling abnormal bin.\n")
  stepName = paste0("call",file.suffix)
  if(any(redo==5)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(length(findJobs(reg))==0){
    abCovCallCases.f <- function(samp, files.df, norm.stats.f, bins.f, stitch.dist, lib.loc){
      library(PopSV, lib.loc=lib.loc)
      load(bins.f)
      call.abnormal.cov(files.df=files.df, samp=samp, out.pdf=paste0(samp,"-sdest-abCovCall.pdf"), FDR.th=.001, merge.cons.bins="stitch", z.th="sdest", norm.stats=norm.stats.f, stitch.dist=stitch.dist, gc.df=bins.df,  min.normal.prop=.6)
    }
    batchMap(reg, abCovCallCases.f, files.df$sample, more.args=list(files.df=files.df, norm.stats.f=out.files[2], bins.f=bins.f, stitch.dist=5e3, lib.loc=lib.loc))
    submitJobs(reg, findJobs(reg))
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  if(status) showStatus(reg)

  if(loose){
    message("\n== 6) Calling abnormal bin with loose threshold.\n")
    stepName = paste0("callLoose",file.suffix)
    if(any(redo==6)) unlink(paste0(stepName, "-files"), recursive=TRUE)
    reg <- makeRegistry(id=stepName, seed=123)
    if(length(findJobs(reg))==0){
      abCovCallCases.f <- function(samp, files.df, norm.stats.f, bins.f, stitch.dist, lib.loc){
        library(PopSV, lib.loc=lib.loc)
        load(bins.f)
        project = subset(files.df, sample==samp)$project
        call.abnormal.cov(files.df=files.df, samp=samp, out.pdf=paste0(project,"/",samp,"/",samp,"-sdest-abCovCall.pdf"), FDR.th=.05, merge.cons.bins="stitch", z.th="sdest", norm.stats=norm.stats.f, stitch.dist=stitch.dist, gc.df=bins.df,  min.normal.prop=.6)
      }
      batchMap(reg, abCovCallCases.f, files.df$sample, more.args=list(files.df=files.df, norm.stats.f=out.files[2], bins.f=bins.f, stitch.dist=5e3, lib.loc=lib.loc))
      submitJobs(reg, findJobs(reg))
      waitForJobs(reg, sleep=sleep)
    }
    if(length(findJobs(reg))!=length(findDone(reg))){
      showStatus(reg)
      if(length(findExpired(reg))>0){
        message("Re-submitting ", findExpired(reg))
        submitJobs(reg, findExpired(reg))
      }
      if(length(findNotSubmitted(reg))>0){
        message("Re-submitting ", findNotSubmitted(reg))
        submitJobs(reg, findNotSubmitted(reg))
      }
      waitForJobs(reg, sleep=sleep)
      if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
    }
    if(status) showStatus(reg)
  }

  res.df = do.call(rbind, reduceResultsList(reg))
  return(res.df)
}
