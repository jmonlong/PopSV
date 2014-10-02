##' Detect abnormal bin from the Z-score distribution. A normal distribution is first
##' fitted to the Z-score distribution. P-values are computed from this estimated
##' null distribution and corrected for multiple testing. Eventually consecutive bins
##' with abnormal read counts can be merged.
##'
##' DETAILS ON THE MERGING APPROACH
##'
##' DETAILS ON THE NORMAL FITTING
##' @title Call abnormal bins
##' @param z the name of the file with the Z-scores for all samples OR or a
##' data.frame with the Z-scores for all samples. 
##' @param sample the name of the sample to analyze.
##' @param out.pdf the name of the output pdf file.
##' @param FDR.th the False Discovery Rate to use for the calls. Further filtering
##' can always be performed.
##' @param merge.cons.bins Should consecutive abnormal bins be merged. If FALSE (default)
##' single bin calls are reported. If TRUE, reported bins are merged. 
##' @param cn the name of the file with the copy number estimates for all samples OR
##' or a data.frame with the copy number estimates for all samples. 
##' @return a data.frame with columns
##' \item{chr, start, end}{the genomic region definition}
##' \item{z}{the Z-score}
##' \item{pv, qv}{the P-value and Q-value(~FDR)}
##' \item{cn.coeff}{the copy number estimate (if 'cn' was not NULL).}
##' \item{nbBinCons}{the number of consecutive bins (if the bins were merged;
##' 'merge.cons.bins=TRUE')}
##' @author Jean Monlong
##' @export
call.abnormal.cov <- function(z,sample,out.pdf=NULL,FDR.th=.05, merge.cons.bins=FALSE, cn=NULL){

    ## load Z-scores and CN coefficients
    if(is.character(z) & length(z)==1){
        headers = read.table(z,nrows=1,as.is=TRUE)
        colC = rep("NULL",length(headers))
        if(!all(c("chr","start","end",sample) %in% headers)){
            stop("Columns missing in Z file. Check that 'chr', 'start', 'end' and the sample column are present.")
        }
        colC[headers %in% c("chr","start","end",sample)] = c("character", rep("integer",2),"numeric")
        res.df = read.table(z,header=TRUE,colClasses=colC)
    } else {
        res.df = z[,c("chr","start","end",sample)]
        rm(z)
    }
    colnames(res.df)[4] = "z"
    if(!is.null(cn)){
        if(is.character(cn) & length(cn)==1){
            headers = read.table(cn,nrows=1,as.is=TRUE)
            colC = rep("NULL",length(headers))
            if(all(headers!=sample)){
                stop("Columns missing in CN file. Check that 'chr', 'start', 'end' and the sample column are present.")
            }
            colC[headers==sample] = "numeric"
            cn = read.table(cn,header=TRUE,colClasses=colC)
        } 
        res.df$cn.coeff = cn[,sample]
        rm(cn)
    }
    
    res.df = subset(res.df, !is.na(z) & !is.infinite(z))
    ## Pvalue/Qvalue estimation
    if(all(is.na(res.df$z))) return(NULL)
    fdr = fdrtool.quantile(res.df$z)
    res.df$pv = fdr$pval
    res.df$qv = fdr$qval

    if(!is.null(out.pdf)){
        pdf(out.pdf,13,10)
    }
    if(!is.null(out.pdf) & any(!is.na(res.df$pv))){
        print(ggplot2::ggplot(subset(res.df,abs(z)<10),ggplot2::aes(x=z)) +
              ggplot2::geom_histogram() + 
              ggplot2::xlab("Z-score") + 
              ggplot2::ylab("number of bins") + 
              ggplot2::theme_bw())
        print(ggplot2::ggplot(res.df,ggplot2::aes(x=pv)) + ggplot2::geom_histogram() +
              ggplot2::xlab("P-value") + 
              ggplot2::ylab("number of bins") + 
              ggplot2::theme_bw())
    }
    
    if(merge.cons.bins){
        res.df = mergeConsBin.z(res.df, fdr.th=FDR.th, sd.null=fdr$sigma.est)

        print(ggplot2::ggplot(res.df,ggplot2::aes(x=factor(nbBinCons))) +
              ggplot2::geom_histogram() +
              ggplot2::ylab("number of bins") + 
              ggplot2::xlab("number of consecutive abnormal bins") +
              ggplot2::theme_bw())

        if(any(colnames(res.df)=="cn.coeff") & any(res.df$nbBinCons>2)) {
            print(ggplot2::ggplot(subset(res.df, nbBinCons>2),ggplot2::aes(x=2*cn.coeff)) +
                  ggplot2::geom_histogram() + ggplot2::theme_bw() +
                  ggplot2::ylab("number of bins") + 
                  ggplot2::xlab("copy number estimate") +
                  ggplot2::ggtitle("At least 3 consecutive abnormal bins") + 
                  ggplot2::xlim(0,5))
        }
    }

    if(!is.null(out.pdf)){
        dev.off()
    }

    if(nrow(res.df)>0 & merge.cons.bins){
        return(data.frame(sample=sample,res.mz))
    } else if(any(res.df$qv <= FDR.th, na.rm=TRUE)){
        return(data.frame(sample=sample,subset(res.df, qv<=FDR.th)))
    } else {
        return(NULL)
    }
}
