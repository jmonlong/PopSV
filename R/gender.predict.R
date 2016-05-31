##' Predict the gender of a sample or a group of samples
##' @title Predict the gender of sample(s)
##' @param bc.df a data.frame with bin coverage. Columns must be "chr", "start", "end". If only one sample is inputed the fourth and last column must be "bc"; if multiple samples are inputed each additionnal column represent one sample.
##' @return a vector or data.frame with the gender for the sample(s)
##' @author Jean Monlong
##' @export
gender.predict <- function(bc.df){
  if(ncol(bc.df)==4){
    bc.df = bc.df[which(bc.df$bc>0),]
    bc.chr = aggregate(bc~chr, data=bc.df, median, na.rm=TRUE)
    bc.chr$bc = bc.chr$bc / median(bc.chr$bc)
    if(bc.chr$bc[which(bc.chr$chr == "X")] < .75){
      return("male")
    } else {
      return("female")
    }
  } else {
    samples = setdiff(colnames(bc.df), c("chr","start","end"))
    covered = apply(bc.df[, samples], 1, function(bc)any(bc>10))
    bc.df = bc.df[which(covered),]
    normBc <- function(bc, chrs){
      bc.chr = tapply(bc, chrs, median, na.rm=TRUE)
      bc.chr = bc.chr / median(bc.chr)
      bc.chr["X"]
    }
    bc.x = apply(bc.df[,samples], normBc, chrs=bc.df$chr)
    bc.x.cl = cutree(hclust(bc.x), 2)
    bc.x.cl.v = tapply(bc.x, bc.x.cl, median)
    return(data.frame(sample=samples, gender=ifelse(bc.x.cl==which.min(bc.x.cl.v),"male","female")))
  }
}
