#' Detect abnormal read coverage using a set of samples as references. From high-throughput sequencing, the genome is binned and read mapping to each bin are counted. From these read counts, PopSV test, in each bin, is the read counts are different between a sample of interest and a set of reference samples. Functions for pre/post-processing are also provided.
#'
#' \tabular{ll}{
#' Package: \tab PopSV\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2015-12-09\cr
#' License: \tab GPL-2\cr
#' }
#' @docType package
#' @name PopSV-package
#' @title Population-based detection of structural variants from Read-Depth signal.
#' @author Jean Monlong \email{jean.monlong@@mail.mcgill.ca}
#' @seealso \url{http://jmonlong.github.io/PopSV/}
#' @useDynLib PopSV
#' @importFrom Rcpp sourceCpp
NULL
