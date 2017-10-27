####################################################################################################################################
### Filename:    f1_sub0.R
### Description: Function for calculating the test statistic for only one whole-plot factor
###              
###
###
####################################################################################################################################

#' Test for main group effect (weighted/unweighted)
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.1.none <- function(X, alpha , group, subject, data, formula ){
  
  temp0 <- hrm.1w.0f(X, alpha , group,  subject, data, "A", paste(as.character(group), " weighted"))
  temp1 <- hrm.1w.0f(X, alpha , group,  subject, data, "Au", paste(as.character(group), " unweighted"))
  
  output <- list()
  output$result <- rbind(temp0, temp1)
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group), c("none"))
  class(output) <- "HRM"
  
  return (output)
}


#' Test for interaction of factor A and B
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the response variable
#' @param H string specifying the hypothesis
#' @param text a string, which will be printed in the output
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.1w.0f <- function(X, alpha, group, subject, data, H, text ){
  
  stopifnot(is.data.frame(X),is.character(subject), is.character(group),alpha<=1, alpha>=0)
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0  
  
  
  group <- as.character(group)
  #factor1 <- as.character(factor1)
  subject <- as.character(subject)
  X <- split(X, X[,group], drop=TRUE)
  a <- length(X)
  d <- 1 #nlevels(X[[1]][,factor1])
  c <- 1
  n <- rep(0,a) 
  
  for(i in 1:a){
    X[[i]] <- X[[i]][ order(X[[i]][,subject]), ]
    X[[i]]<-X[[i]][,data]
    X[[i]] <- matrix(X[[i]],ncol=d*c,byrow=TRUE)
    n[i] <- length(X[[i]])
  }
  
  # creating X_bar (list with a entries)
  X_bar <- sapply(X, colMeans, na.rm=TRUE)
  
  
  if(H=="A"){
    K <- 1
    S <- diag(n)-1/sum(n)*n%*%t(n)
  } else if(H=="Au"){
    K <- 1
    S <- P(a)
  }
  
  
  # creating dual empirical covariance matrices
  K_A <- kronecker(S, K)
  V <- lapply(X, DualEmpirical2, B=K)
  
  #################################################################################################
  
  # f
  f_1 <- 0
  f_2 <- 0
  
  for(i in 1:a){
    f_1 <- f_1 + (S[i,i]*1/n[i])^2*.E1(n,i,V[[i]])
    j <- i+1
    while(j<=a){
      f_1 <- f_1 + 2*(S[i,i]*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j <- j+1
    }
  }
  
  for(i in 1:a){
    f_2 <- f_2 + (S[i,i]*1/n[i])^2*.E2(n,i,V[[i]])
    j <- i+1
    while(j<=a){
      f_2 <- f_2 + 2*S[i,j]*S[j,i]*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*K%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%K%*%t(X[[i]])%*%P(n[i]))
      j <- j+1
    }
  }
  
  f <- f_1/f_2
  
  
  ##################################################################################################
  
  
  
  #################################################################################################
  # f0
  f0_1 <- f_1
  f0_2 <- 0
  
  
  for(i in 1:a){
    f0_2 <- f0_2 + (S[i,i]*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0 <- f0_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit <- qf(1-alpha,f,f0)
  
  # Test
  direct <- direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct <- direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  
  test <- (t(X_bar)%*%K_A%*%X_bar)/(t(rep(1,dim(K_A)[1]))%*%(K_A*direct)%*%(rep(1,dim(K_A)[1])))
  p.value <- 1-pf(test,f,f0)
  output <- data.frame(hypothesis=text,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# Hypothesis AC End ------------------------------------------------------------
