##' When two cohort have been run separately but the variants have to be compared, this function finds the bins in which both runs had comparable detection power. Or at least it removes most of the bins that are clearly not covered similarly in both runs. H
##'
##' If 'plot=TRUE' (Default), a graph is displayed. It compares the average coverage in both cohort and shows the thresholds used to filter out bins with different coverage.
##' @title Find comparable bins between two sample cohort.
##' @return a data.frame with the bins to use for the comparison.
##' @author Jean Monlong
##' @export
##' @param msd1 a data.frame or path to the file with the normalization stats from the first run.
##' @param msd2 a data.frame or path to the file with the normalization stats from the second run.
##' @param plot should a graph be displayed. Default is TRUE.
find.comparable.bins <- function(msd1, msd2, plot=TRUE){
  ## Uglily appease R checks
  chr = start = NULL

  ## Read the files if necessary
  if(!is.data.frame(msd1)){
    if(length(msd1)>1){
      stop("'msd1' should be either a data.frame or a character vector of length 1 with the path to the correct file.")
    }
    if(!file.exists(msd1)){
      stop("msd1: ",msd1, " file does not exist.")
    }
    msd1 = data.table::fread(msd1)
    data.table::setkey(msd1, chr, start)
    msd1 = as.data.frame(msd1)
  }
  if(!is.data.frame(msd2)){
    if(length(msd2)>1){
      stop("'msd2' should be either a data.frame or a character vector of length 1 with the path to the correct file.")
    }
    if(!file.exists(msd2)){
      stop("msd2: ",msd2, " file does not exist.")
    }
    msd2 = data.table::fread(msd2)
    data.table::setkey(msd2, chr, start)
    msd2 = as.data.frame(msd2)
  }

  ## Merge them
  cols = c("chr","start","end","m")
  m.df = merge(msd1[,cols], msd2[,cols], by=cols[1:3], suffixes = 1:2)
  m.df = m.df[which(!is.na(m.df$m1) & !is.na(m.df$m2)),]

  ## Polar Kmeans
  theta = atan(log10(m.df$m2+1)/log10(m.df$m1+1))
  km = kmeans(theta, 4)
  extreme.cl = c(which.max(km$centers), which.min(km$centers))
  m.df$comparable = !(km$cluster %in% extreme.cl)

  ## Graph
  if(plot){
    ## Uglily appease R checks
    m1 = m2 = ..count.. = NULL
    th = sort(km$centers)
    th = c(mean(c(th[1],th[2])), mean(c(th[3],th[4])))
    bw = c(mean(log10(1+m.df$m1)), mean(log10(1+m.df$m2))) /100
    print(ggplot2::ggplot(m.df, ggplot2::aes(x=m1+1,y=m2+1, fill=log10(..count..))) +
          ggplot2::geom_bin2d(binwidth=bw) + ggplot2::scale_x_log10() +
          ggplot2::scale_y_log10() + ggplot2::geom_abline(slope=tan(th), linetype=2) +
          ggplot2::theme_bw() + ggplot2::xlab("average coverage 1") +
          ggplot2::ylab("average coverage 2") +
          ggplot2::scale_fill_gradient(name="nb bins\n(log10)"))
  }

  ## Return comparable bins
  m.df = m.df[which(m.df$comparable),]
  m.df$m1 = m.df$m2 = m.df$comparable = NULL
  return(m.df)  
}
