call.abnormal.cov.multiSamp <- function(z,samples,out.pdf=NULL,FDR.th=.05, chunk.size=100000){
    ## load Z-scores and FC coefficients
    if(is.character(z) & length(z)==1){
        res.l = list()
        con = file(z,"r")
        headers = unlist(strsplit(readLines(con,n=1),"\t"))
        if(!all(c("chr","start","end",samples) %in% headers)){
            stop("Columns missing in Z file. Check that 'chr', 'start', 'end' and the sample columns are present.")
        }
        while(length((lines = readLines(con,n=chunk.size)))>0){
            z.chunk = matrix(unlist(strsplit(lines, "\t")), ncol=length(headers), byrow=TRUE)
            res.chunk = data.frame(chr=z.chunk[,1], start=as.integer(z.chunk[,2]), end=as.integer(z.chunk[,3]), stringsAsFactors=FALSE)
            z.chunk = apply(z.chunk[,headers%in%samples], 2, function(ee) norm.z.quantile(as.numeric(ee)))
            res.chunk$chi=sqrt(rowSums(z.chunk*z.chunk))
            res.l = c(res.l, list(res.chunk))
        }
        res.df = plyr::ldply(res.l, identity)
        rm(res.l)
        close(con)
    } else {
        res.df = data.frame(z[,c("chr","start","end")], chi=sqrt(rowSums(z[,samples]*z[,samples])), stringsAsFactors=FALSE)
        rm(z)
    }
    
    res.df = res.df[which(!is.na(res.df$chi)),]
    ## Pvalue/Qvalue estimation
    if(all(is.na(res.df$chi))) return(NULL)
    res.df$pv = pchisq(res.df$chi^2, df=length(samples), lower.tail=FALSE)
    if(any(res.df$pv<1e-50)) res.df$pv[res.df$pv==0] = 1e-50
    ft = fdrtool::fdrtool(res.df$pv,statistic="pvalue",plot=FALSE,verbose=FALSE)
    if(any(ft$qval<.05,na.rm=TRUE) | !any(res.df$pv<1e-10)){
        res.df$qv = ft$qval
    } else {
        res.df$qv = p.adjust(res.df$pv,method="fdr")
    }
    
    if(!is.null(out.pdf)){
        pdf(out.pdf,13,10)
    }
    if(!is.null(out.pdf) & any(!is.na(res.df$pv))){
      pv = NULL ## Uglily appease R checks
        print(ggplot2::ggplot(res.df,ggplot2::aes(x=pv)) + ggplot2::geom_histogram() +
              ggplot2::xlab("P-value") + ggplot2::xlim(0,1) + 
              ggplot2::ylab("number of bins") + 
              ggplot2::theme_bw())
        print(ggplot2::ggplot(res.df,ggplot2::aes(x=ppoints(length(pv)),y=sort(pv))) + ggplot2::geom_point() + ggplot2::scale_y_log10() + ggplot2::scale_x_log10() +
              ggplot2::xlab("Expected P-value") + 
              ggplot2::ylab("Observed P-value") + 
              ggplot2::theme_bw())
    }
    
    if(!is.null(out.pdf)){
        dev.off()
    }

}
