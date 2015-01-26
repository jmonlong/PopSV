##' Winsorising a vector: any value higher(lower) than U(L) are changed to U(L).
##' @title Winsorize a vector
##' @param x the input to winsorize.
##' @param u the highest value desired.
##' @param l the lowest value desired.
##' @return a vector with winsorized values.
##' @author Jean Monlong
##' @keywords internal
winsor <- function(x, u=5, l=NULL){
  if(!is.null(u) & any(x>u)){
    x[x>u] = u
  }
  if(!is.null(l) & any(x>l)){
    x[x<l] = l
  }
  return(x)
}
