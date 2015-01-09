breakpoint.finder <- function(bkpt.gr, files.df, test.sample, ref.samples, bp.res=10, proper=TRUE, map.quality=30, slidW.size=300,slidW.step=2){
    bkpt.gr = GenomicRanges::split(bkpt.gr, 1:length(bkpt.gr))

    gr.bam <- function(bam.file,gr,proper=TRUE,map.quality=30){
        bai.file = sub("bam$","bai",bam.file,perl=TRUE)
        if(!file.exists(bai.file)){
            bai.file = paste0(bam.file,".bai")
            if(!file.exists(bai.file)){
                stop("Index file is missing (neither '.bai' nor '.bam.bai').")
            }
        }
        param = Rsamtools::ScanBamParam(which=gr,
            what = c("rname", "pos", "qwidth","mapq"),
            flag = Rsamtools::scanBamFlag(isProperPair=proper,isDuplicate=FALSE,isNotPassingQualityControls=FALSE,isUnmappedQuery=FALSE)
            )
        bam = Rsamtools::scanBam(bam.file, index=bai.file,param=param)
        bam.df = plyr::ldply(bam, as.data.frame, .id="bin")
        bam.df = subset(bam.df, mapq>= map.quality)
        return(with(bam.df, GenomicRanges::GRanges(rname, IRanges::IRanges(pos, width=qwidth), bin=bin)))
    }

    cov.gr <- function(gr.i,gr.i.frag){
        gr.cov = GRanges(seqnames(gr.i)[1], IRanges(gr.bk[-length(gr.bk)], width=bin.size))
        cov.t = table(cut((start(gr.i)+end(gr.i))/2, gr.bk, labels=FALSE))
        gr.cov$cov = 0
        gr.cov$cov[as.integer(names(cov.t))] = as.integer(cov.t)
        gr.cov
    }

    comp.diff <- function(cov.m){
        cov.m = apply(cov.m, 2, function(ee)ee/median(ee))
        median(apply(cov.m[,-1], 2, function(ee)sum(abs(ee-cov.m[,1]))))
    }
    
    find.bp <- function(gr.f){
        if(bp.res>1){
            gr.bk = seq(start(gr.f), end(gr.f), bp.res)
            gr.frag = GRanges(seqnames(gr.f)[1], IRanges(gr.bk[-length(gr.bk)], width=bp.res))
            cov.l = lapply(c(test.sample, ref.samples), function(samp.i){
                gr.i = gr.bam(files.df$bam[files.df$sample==samp.i], gr.f, proper=proper, map.quality=map.quality)
                cov.v = rep(0, length(gr.frag))
                cov.t = table(cut((start(gr.i)+end(gr.i))/2, gr.bk, labels=FALSE))
                cov.v[as.integer(names(cov.t))] = as.integer(cov.t)
                cov.v
            })
            cov = matrix(0, length(gr.frag), length(ref.samples)+1)
            for(ii in 1:(length(ref.samples)+1)) cov[,ii] = cov.l[[ii]]
            diff.v = sapply(seq(1,length(gr.frag)-slidW.size+1,slidW.step), function(sw.ii){
                comp.diff(cov[sw.ii:(sw.ii+slidW.size-1),])
            })
            return(gr.frag[which.max(diff.v)*slidW.step+slidW.size/2])
        } else {
            stop("Single nucleotide resolution not implemented yet...")
        }
    }

    bkpt.gr = lapply(bkpt.gr, find.bp)

    return(bkpt.gr)
}
