##' Create an empty data.frame with NA and specified column types.
##' @title Create empty data.frame
##' @param colClasses a vector with the types definition for the columns.
##' @param nrows the number of rows in the data.frame to create.
##' @return a data.frame with 'NA' but of the specified types
##' @author Jean Monlong
##' @keywords internal
createEmptyDF <- function(colClasses, nrows){
    df = as.data.frame(matrix(NA, nrows,length(colClasses)))
    for(ii in 1:length(colClasses)) class(df[,ii]) = colClasses[ii]
    df
}
