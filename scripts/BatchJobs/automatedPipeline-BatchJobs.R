message(" /!\ Might need to be tweaked /!\ ")

message("Functions :
- 'autoGCcounts' to count BC in each sample.
- 'autoNormTest' to normalize and test all the samples.
- 'autoExtra' for some other functions.
")

autoGCcounts <- function(files.f,
                         bins.f,
                         redo=NULL,
                         sleep=180,
                         status=FALSE,
                         file.suffix="",
                         lib.loc=NULL,
                         other.resources=NULL,
                         skip=NULL,
                         step.walltime=c(2,20),
                         step.cores=c(1,1),
                         resetExpired=FALSE){
  load(files.f)
  step.walltime = paste0(step.walltime, ":0:0")

  message("\n== 1) Get GC content in each bin.\n")
  stepName = paste0("getGC",file.suffix)
  if(any(redo==1)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(!any(skip==1) & length(findJobs(reg))==0){
    getGC.f <- function(imF){
      load(imF)
      library(PopSV, lib.loc=lib.loc)
      bins.df = getGC.hg19(bins.df)
      save(bins.df, file=imF)
    }
    batchMap(reg, getGC.f,bins.f)
    submitJobs(reg, findJobs(reg), resources=c(list(walltime=step.walltime[1], nodes="1", cores=step.cores[1]), other.resources))
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[1], nodes="1", cores=step.cores[1]), other.resources))
    }
    if(resetExpired){
      resetJobs(reg, c(findErrors(reg), findExpired(reg)))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[1], nodes="1", cores=step.cores[1]), other.resources))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  if(status) showStatus(reg)

  message("\n== 2) Get bin counts in each sample and correct for GC bias.\n")
  stepName = paste0("getBC",file.suffix)
  if(any(redo==2)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(!any(skip==2) & length(findJobs(reg))==0){
    getBC.f <- function(file.i, bins.f, files.df, lib.loc){
      library(PopSV, lib.loc=lib.loc)
      load(bins.f)
      bam.f = files.df$bam[file.i]
      if("bam2" %in% colnames(files.df)) bam.f = c(bam.f, files.df$bam2[file.i])
      bb.o = bin.bam(bam.f, bins.df, files.df$bc[file.i])
      correct.GC(files.df$bc.gz[file.i], bins.df, files.df$bc.gc[file.i])
      bb.o
    }
    batchMap(reg, getBC.f,1:nrow(files.df), more.args=list(bins.f=bins.f, files.df=files.df, lib.loc=lib.loc))
    submitJobs(reg, findJobs(reg), resources=c(list(walltime=step.walltime[2], nodes="1", cores=step.cores[2]), other.resources))
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[2], nodes="1", cores=step.cores[2]), other.resources))
    }
    if(resetExpired){
      resetJobs(reg, c(findErrors(reg), findExpired(reg)))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[2], nodes="1", cores=step.cores[2]), other.resources))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  if(status) showStatus(reg)
}

autoNormTest <- function(files.f,
                         bins.f,
                         redo=NULL,
                         rewrite=FALSE,
                         sleep=180,
                         status=FALSE,
                         loose=FALSE,
                         file.suffix="",
                         lib.loc=NULL,
                         other.resources=NULL,
                         norm=c("1pass","trim"),
                         ref.samples=NULL,
                         FDR.th=.001,
                         col.bc="bc.gc.gz",
                         step.walltime=c(10,12,6,10,3,3),
                         step.cores=c(12,2,3,1,1,1),
                         file.suffix.ref=NULL,
                         skip=NULL,
                         resetExpired=FALSE){
  if(is.null(file.suffix.ref)){
    file.suffix.ref = file.suffix
  }
  load(files.f)
  step.walltime = paste0(step.walltime, ":0:0")

  message("\n== 1) Sample QC and reference definition.\n")
  bc.ref.f = paste0("bc-gcCor",file.suffix.ref,".tsv")
  sampQC.pdf.f = paste0("sampQC",file.suffix.ref,".pdf")
  stepName = paste0("sampQC",file.suffix.ref)
  if(any(redo==1)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(!any(skip==1) & length(findJobs(reg))==0){
    if(!is.null(ref.samples)){
      files.ref = subset(files.df, sample %in% ref.samples)
    } else {
      files.ref = files.df
    }
    sampQC.f <- function(bc.all.f, bins.f, files.df, sampQC.pdf.f, lib.loc, col.bc){
      load(bins.f)
      library(PopSV, lib.loc=lib.loc)
      pdf(sampQC.pdf.f)
      qc.o = qc.samples(files.df, bins.df, bc.all.f, nb.cores=6, nb.ref.samples=200, col.bc=col.bc)
      dev.off()
      qc.o
    }
    batchMap(reg, sampQC.f,bc.ref.f, more.args=list(bins.f=bins.f, files.df=files.ref, sampQC.pdf.f=sampQC.pdf.f, lib.loc=lib.loc, col.bc=col.bc))
    submitJobs(reg, 1, resources=c(list(walltime=step.walltime[1], nodes="1", cores=step.cores[1]), other.resources))
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[1], nodes="1", cores=step.cores[1]), other.resources))
    }
    if(resetExpired){
      resetJobs(reg, c(findErrors(reg), findExpired(reg)))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[1], nodes="1", cores=step.cores[1]), other.resources))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  samp.qc.o = loadResult(reg, 1)
  save(samp.qc.o, file=paste0(stepName,".RData"))
  if(status) showStatus(reg)

  message("\n== 2) Reference sample normalization.\n")
  stepName = paste0("bcNormTN",file.suffix.ref)
  if(any(redo==2)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(!any(skip==2) & length(findJobs(reg))==0){
    load(bins.f)
    if(all(colnames(bins.df)!="sm.chunk")){
      bins.df = chunk.bin(bins.df, bg.chunk.size=5e5, sm.chunk.size=1e4)
      save(bins.df, file=bins.f)
    }
    bcNormTN.f <- function(chunk.id, file.bc, file.bin, cont.sample, lib.loc, norm){
      load(file.bin)
      library(PopSV, lib.loc=lib.loc)
      bc.df = read.bedix(file.bc, subset(bins.df, bg.chunk==subset(bins.df, sm.chunk==chunk.id)$bg.chunk[1]))
      tn.norm(bc.df, cont.sample, bins=subset(bins.df, sm.chunk==chunk.id)$bin, norm=norm, force.diff.chr=TRUE)
    }
    batchMap(reg, bcNormTN.f,unique(bins.df$sm.chunk), more.args=list(file.bc=samp.qc.o$bc, file.bin=bins.f,cont.sample=samp.qc.o$cont.sample, lib.loc=lib.loc, norm=norm))
    submitJobs(reg, findJobs(reg) , resources=c(list(walltime=step.walltime[2], nodes="1", cores=step.cores[2]), other.resources))
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[2], nodes="1", cores=step.cores[2]), other.resources))
    }
    if(resetExpired){
      resetJobs(reg, c(findErrors(reg), findExpired(reg)))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[2], nodes="1", cores=step.cores[2]), other.resources))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  ## Write normalized bin counts and reference metrics
  out.files = paste(paste0("ref",file.suffix.ref), c("bc-norm.tsv", "norm-stats.tsv"), sep="-")
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
  stepName = paste0("zRef",file.suffix.ref)
  if(any(redo==3)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(!any(skip==3) & length(findJobs(reg))==0){
    files.ref = subset(files.df, sample %in% samp.qc.o$ref.samples)
    files.l = tapply(1:nrow(files.ref), cut(1:nrow(files.ref), nrow(files.ref)/5), function(ii)files.ref[ii,])
    zRef.f <- function(files.ii, bc.f, files.l, ns.f, lib.loc, nb.cores){
      library(PopSV, lib.loc=lib.loc)
      z.comp(bc.f=bc.f, norm.stats.f=ns.f, files.df=files.l[[files.ii]], nb.cores=nb.cores, z.poisson=TRUE, chunk.size=1e4)
    }
    batchMap(reg, zRef.f,1:length(files.l), more.args=list(bc.f=out.files[1], files.l=files.l, ns.f=out.files[2], lib.loc=lib.loc, nb.cores=step.cores[3]))
    submitJobs(reg, findJobs(reg), resources=c(list(walltime=step.walltime[3], nodes="1", cores=step.cores[3]), other.resources))
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[3], nodes="1", cores=step.cores[3]), other.resources))
    }
    if(resetExpired){
      resetJobs(reg, c(findErrors(reg), findExpired(reg)))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[3], nodes="1", cores=step.cores[3]), other.resources))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  if(status) showStatus(reg)

  message("\n== 4) Normalization and Z-score computation for other samples.\n")
  stepName = paste0("zOthers",file.suffix)
  if(any(redo==4)) unlink(paste0(stepName, "-files"), recursive=TRUE)
  reg <- makeRegistry(id=stepName, seed=123)
  if(!any(skip==4) & length(findJobs(reg))==0){
    callOthers.f <- function(samp, cont.sample, files.df, norm.stats.f, bc.ref.f, lib.loc, col.bc){
      library(PopSV, lib.loc=lib.loc)
      tn.test.sample(samp, files.df, cont.sample, bc.ref.f, norm.stats.f, z.poisson=TRUE, aberrant.cases=FALSE, col.file=col.bc)
    }
    batchMap(reg, callOthers.f,setdiff(files.df$sample, samp.qc.o$ref.samples), more.args=list(cont.sample=samp.qc.o$cont.sample, files.df=files.df, norm.stats.f=out.files[2], bc.ref.f=samp.qc.o$bc, lib.loc=lib.loc, col.bc=col.bc))
    submitJobs(reg, findJobs(reg), resources=c(list(walltime=step.walltime[4], nodes="1", cores=step.cores[4]), other.resources))
    waitForJobs(reg, sleep=sleep)
  }
  if(length(findJobs(reg))!=length(findDone(reg))){
    showStatus(reg)
    if(length(findExpired(reg))>0){
      message("Re-submitting ", findExpired(reg))
      submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[4], nodes="1", cores=step.cores[4]), other.resources))
    }
    if(resetExpired){
      resetJobs(reg, c(findErrors(reg), findExpired(reg)))
    }
    if(length(findNotSubmitted(reg))>0){
      message("Re-submitting ", findNotSubmitted(reg))
      submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[4], nodes="1", cores=step.cores[4]), other.resources))
    }
    waitForJobs(reg, sleep=sleep)
    if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
  }
  if(status) showStatus(reg)

  if(!loose){
    message("\n== 5) Calling abnormal bin.\n")
    stepName = paste0("call",file.suffix)
    if(any(redo==5)) unlink(paste0(stepName, "-files"), recursive=TRUE)
    reg <- makeRegistry(id=stepName, seed=123)
    if(!any(skip==5) & length(findJobs(reg))==0){
      abCovCallCases.f <- function(samp, files.df, norm.stats.f, bins.f, stitch.dist, lib.loc, FDR.th){
        library(PopSV, lib.loc=lib.loc)
        load(bins.f)
        call.abnormal.cov(files.df=files.df, samp=samp, out.pdf=paste0(samp,"-sdest-abCovCall.pdf"), FDR.th=FDR.th, merge.cons.bins="cbs", z.th="sdest", norm.stats=norm.stats.f, stitch.dist=stitch.dist, gc.df=bins.df,  min.normal.prop=.6, sub.z=1e3)
      }
      batchMap(reg, abCovCallCases.f, files.df$sample, more.args=list(files.df=files.df, norm.stats.f=out.files[2], bins.f=bins.f, stitch.dist=5e3, lib.loc=lib.loc, FDR.th=FDR.th))
      submitJobs(reg, findJobs(reg) , resources=c(list(walltime=step.walltime[5], nodes="1", cores=step.cores[5]), other.resources))
      waitForJobs(reg, sleep=sleep)
    }
    if(length(findJobs(reg))!=length(findDone(reg))){
      showStatus(reg)
      if(length(findExpired(reg))>0){
        message("Re-submitting ", findExpired(reg))
        submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[5], nodes="1", cores=step.cores[5]), other.resources))
      }
      if(resetExpired){
        resetJobs(reg, c(findErrors(reg), findExpired(reg)))
      }
      if(length(findNotSubmitted(reg))>0){
        message("Re-submitting ", findNotSubmitted(reg))
        submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[5], nodes="1", cores=step.cores[5]), other.resources))
      }
      waitForJobs(reg, sleep=sleep)
      if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
    }
    if(status) showStatus(reg)

  } else {
    message("\n== 6) Calling abnormal bin with loose threshold.\n")
    stepName = paste0("callLoose",file.suffix)
    if(any(redo==6)) unlink(paste0(stepName, "-files"), recursive=TRUE)
    reg <- makeRegistry(id=stepName, seed=123)
    if(!any(skip==6) & length(findJobs(reg))==0){
      abCovCallCases.f <- function(samp, files.df, norm.stats.f, bins.f, stitch.dist, lib.loc){
        library(PopSV, lib.loc=lib.loc)
        load(bins.f)
        project = subset(files.df, sample==samp)$project
        call.abnormal.cov(files.df=files.df, samp=samp, out.pdf=paste0(samp,"/",samp,"-sdest-abCovCall.pdf"), FDR.th=.05, merge.cons.bins="stitch", z.th="sdest", norm.stats=norm.stats.f, stitch.dist=stitch.dist, gc.df=bins.df,  min.normal.prop=.6)
      }
      batchMap(reg, abCovCallCases.f, files.df$sample, more.args=list(files.df=files.df, norm.stats.f=out.files[2], bins.f=bins.f, stitch.dist=5e3, lib.loc=lib.loc))
      submitJobs(reg, findJobs(reg) , resources=c(list(walltime=step.walltime[6], nodes="1", cores=step.cores[6]), other.resources))
      waitForJobs(reg, sleep=sleep)
    }
    if(length(findJobs(reg))!=length(findDone(reg))){
      showStatus(reg)
      if(length(findExpired(reg))>0){
        message("Re-submitting ", findExpired(reg))
        submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[6], nodes="1", cores=step.cores[6]), other.resources))
      }
      if(resetExpired){
        resetJobs(reg, c(findErrors(reg), findExpired(reg)))
      }
      if(length(findNotSubmitted(reg))>0){
        message("Re-submitting ", findNotSubmitted(reg))
        submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[6], nodes="1", cores=step.cores[6]), other.resources))
      }
      waitForJobs(reg, sleep=sleep)
      if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
    }
    if(status) showStatus(reg)
  }

  res.df = do.call(rbind, reduceResultsList(reg))
  return(res.df)
}


autoExtra <- function(files.f,
                      bins.f,
                      do=NULL,
                      redo=NULL,
                      sleep=180,
                      status=FALSE,
                      file.suffix="",
                      lib.loc=NULL,
                      other.resources=NULL,
                      step.walltime=c(6,6),
                      step.cores=c(3,2),
                      col.files="bc.gc.gz",
                      bc.ref.f="ref-bc-norm.tsv",
                      resetExpired=FALSE,
                      seed.c=123){
  load(files.f)
  step.walltime = paste0(step.walltime, ":0:0")

  if(is.null(do)){
    stop("Option for 'do=': '1' quick counts, '2' split ref bc norm.")
  }

  if(do==1){
    message("\n== 1) Quick counts.\n")
    stepName = paste0("quickCount",file.suffix)
    if(any(redo==1)) unlink(paste0(stepName, "-files"), recursive=TRUE)
    reg <- makeRegistry(id=stepName, seed=seed.c)
    if(length(findJobs(reg))==0){
      quickCount.f <- function(col.files, bins.f, files.df, lib.loc, nb.cores){
        load(bins.f)
        library(PopSV, lib.loc=lib.loc)
        quick.count(files.df, bins.df, col.files=col.files, nb.rand.bins=1e3, nb.cores=nb.cores)
      }
      batchMap(reg, quickCount.f,col.files, more.args=list(bins.f=bins.f, files.df=files.df, lib.loc=lib.loc, nb.cores=step.cores[1]))
      submitJobs(reg, 1, resources=c(list(walltime=step.walltime[1], nodes="1", cores=step.cores[1]), other.resources))
      waitForJobs(reg, sleep=sleep)
    }
    if(length(findJobs(reg))!=length(findDone(reg))){
      showStatus(reg)
      if(length(findExpired(reg))>0){
        message("Re-submitting ", findExpired(reg))
        submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[1], nodes="1", cores=step.cores[1]), other.resources))
      }
      if(resetExpired){
        resetJobs(reg, c(findErrors(reg), findExpired(reg)))
      }
      if(length(findNotSubmitted(reg))>0){
        message("Re-submitting ", findNotSubmitted(reg))
        submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[1], nodes="1", cores=step.cores[1]), other.resources))
      }
      waitForJobs(reg, sleep=sleep)
      if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
    }
    if(status) showStatus(reg)
    return(loadResult(reg, 1))
  }

  if(do==2){
    if(!file.exists(bc.ref.f)){
      stop(bc.ref.f, " not found. Check 'bc.ref.f=' parameter.")
    }
    message("\n== 2) Split BC norm in ref samples.\n")
    stepName = paste0("splitRef",file.suffix)
    if(any(redo==2)) unlink(paste0(stepName, "-files"), recursive=TRUE)
    reg <- makeRegistry(id=stepName, seed=seed.c)
    if(length(findJobs(reg))==0){
      splitRef.f <- function(bc.ref, files.df, lib.loc){
        library(PopSV, lib.loc=lib.loc)
        write.split.samples(list(bc=bc.ref.f), files.df, files.col="bc.gc.norm", reorder=TRUE)
      }
      batchMap(reg, splitRef.f, bc.ref.f, more.args=list(files.df=files.df, lib.loc=lib.loc))
      submitJobs(reg, 1, resources=c(list(walltime=step.walltime[2], nodes="1", cores=step.cores[2]), other.resources))
      waitForJobs(reg, sleep=sleep)
    }
    if(length(findJobs(reg))!=length(findDone(reg))){
      showStatus(reg)
      if(length(findExpired(reg))>0){
        message("Re-submitting ", findExpired(reg))
        submitJobs(reg, findExpired(reg), resources=c(list(walltime=step.walltime[2], nodes="1", cores=step.cores[2]), other.resources))
      }
      if(resetExpired){
        resetJobs(reg, c(findErrors(reg), findExpired(reg)))
      }
      if(length(findNotSubmitted(reg))>0){
        message("Re-submitting ", findNotSubmitted(reg))
        submitJobs(reg, findNotSubmitted(reg), resources=c(list(walltime=step.walltime[2], nodes="1", cores=step.cores[2]), other.resources))
      }
      waitForJobs(reg, sleep=sleep)
      if(length(findJobs(reg))!=length(findDone(reg))) stop("Not done yet or failed, see for yourself")
    }
    if(status) showStatus(reg)
    return(loadResult(reg, 1))
  }
}
