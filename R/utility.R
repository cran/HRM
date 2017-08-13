####################################################################################################################################
### Filename:    utility.R
### Description: trace estimator functions, dual empirical covariance matrix function and print/summary are defined here
###              
###
###
####################################################################################################################################


#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E1 <- function(n,i, M) {
  return ((n[i]*(n[i]-1))/((n[i]-2)*(n[i]+1))*(matrix.trace(M)^2-2/n[i]*matrix.trace(M%*%M)))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E2 <- function(n,i, M) {
  return ((n[i]-1)^2/((n[i]-2)*(n[i]+1))*(matrix.trace(M%*%M)-1/(n[i]-1)*matrix.trace(M)^2))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E3 <- function(M_i, M_j) {
  return (matrix.trace(M_i)*matrix.trace(M_j))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E4 <- function(M_i,M_j) {
  return (matrix.trace(M_i%*%M_j))
}


#' Function for the output: significant p-values have on or more stars
#' 
#' @param value p-value
#' @keywords internal
.hrm.sigcode <- function(value) {
  
  char <- ""
  if(value <= 0.1 & value > 0.05) { char = "."}
  if(value <= 0.05 & value > 0.01) {char = '*'}
  if(value <= 0.01 & value > 0.001) {char = "**"}
  if(value <= 0.001 & value >= 0) {char = "***"}
  return (char)
}

#' @keywords export
print.HRM <- function(x, ...) {
  if(!is.null(x$formula)) {
    cat("Call:", "\n")
    print(x$formula)
    cat("\n")
  }

  print(x$result, row.names = FALSE)
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}
#' @keywords export
summary.HRM <- function(object, ...) {
  cat("Summary:\n")
  cat("\n")
  if(!is.null(object$formula)) {
    cat("Call:", "\n")
    print(object$formula)
    cat("\n")
    cat("between-subject factors: ")
    cat(object$factors[[1]], sep=", ")
    cat("\n")
    cat("within-subject factors: ")
    cat(object$factors[[2]], sep=", ")
    cat("\n")
    cat("\n")
  }
  
  print(object$result, row.names = FALSE)
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}

#' Function for the indentity matrix
#' 
#' @param size dimension of the matrix
#' @keywords internal
I <- function(size){
  return (diag(rep(1,size)))
}
#' Function for a matrix with entries 1
#' 
#' @param size dimension of the matrix
#' @keywords internal

J <- function(size){
  return (rep(1,size)%*%t(rep(1,size)))
}
#' Function for the centering matrix
#' 
#' @param size dimension of the matrix
#' @keywords internal
P <- function(size){
  return (I(size)-J(size)*1/size)
}

#' Function for the dual empirical matrix
#' 
#' @param Data data.frame
#' @param B not used
#' @keywords internal
DualEmpirical <- function(Data, B){
  n <- dim(Data)[1]
  B <- B(dim(Data)[2]) # B = e.g. t(J_d)%*%J_d
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}
#' Function for the dual empirical matrix
#' 
#' @param Data data.frame
#' @param B part of the hypothesis matrix
#' @keywords internal
DualEmpirical2 <- function(Data, B){
  n <- dim(Data)[1]
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}