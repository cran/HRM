####################################################################################################################################
### Filename:    f3_sub3.R
### Description: Functions for calculating the test statistic for only three subplot factor
###              
###
###
####################################################################################################################################

#' Test for two subplot factors
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param factor3 column name of the data frame X of the third factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.3.three <- function(X, alpha , factor1, factor2, factor3, subject, data, formula, testing = rep(1,7), nonparametric ){
  
  ranked <- NULL
  varQGlobal <- NULL
  
  # create list for storing results; NULL used, because it is ignored by rbind
  temp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  for(i in 1:7){
    if(testing[i]) {
      temp[[i]] <-  hrm.0w.3s(X, alpha , factor1,  factor2, factor3, subject, data, i, "", nonparametric, ranked, varQGlobal)
    }
  }
  

  output <- list()
  output$result <- rbind(temp[[1]], temp[[2]], temp[[3]], temp[[4]], temp[[5]], temp[[6]], temp[[7]])
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c("none"), c(factor1, factor2, factor3))
  output$data <- X
  output$var <- varQGlobal
  output$nonparametric <- nonparametric
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
hrm.0w.3s <- function(X, alpha , factor1, factor2, factor3, subject, data, H = 1, text ="", nonparametric, ranked, varQGlobal ){
  
  stopifnot(is.data.frame(X),is.character(subject), is.character(factor1), is.character(factor2), alpha<=1, alpha>=0, is.logical(nonparametric))
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0  
  
  factor1 <- as.character(factor1)
  factor2 <- as.character(factor2)
  factor3 <- as.character(factor3)
  subject <- as.character(subject)
  
  X <- as.data.table(X)
  setnames(X, c(data, factor1, subject, factor2, factor3), c("data", "factor1", "subject", "factor2", "factor3"))
  
  a <- 1
  d <- nlevels(X[,factor1])
  c <- nlevels(X[,factor2])
  c2 <- nlevels(X[,factor3])
  n <- dim(X)[1]
  
  if(nonparametric & is.null(ranked)) {
    X[,data:= 1/(sum(n)*d*c*c2)*(rank(X[,data], ties.method = "average") - 1/2)]
  }
  
  for(i in 1:a){
    X <- X[ order(subject, factor1, factor2, factor3), ]
    X <- X[,data]
    X <- matrix(X,ncol=d*c*c2,byrow=TRUE)
    n[i] <- dim(X)[1]
  }
  
  if(is.null(ranked)){
    eval.parent(substitute(ranked<-X))
  } else {
    X <- ranked
  }
  
  # creating X_bar (list with a entries)
  X_bar <- colMeans(X)  # as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  
  
  
  # main effects
  if(H==1){
    K <- kronecker(kronecker(P(d), 1/c*J(c)), 1/c2*J(c2))
    text <- paste(as.character(factor1) )
  }
  if(H==2){
    K <- kronecker(kronecker(1/d*J(d), P(c)), 1/c2*J(c2))
    text <- paste(as.character(factor2) )
  }
  if(H==3){
    K <- kronecker(kronecker(1/d*J(d), 1/c*J(c)), P(c2))
    text <- paste(as.character(factor3) )
  }
  # inferaction effects of 2 factors
  if(H==4){
    K <- kronecker(kronecker(P(d), P(c)), 1/c2*J(c2))
    text <- paste(as.character(factor1), ":", as.character(factor2) )
  }
  if(H==5){
    K <- kronecker(kronecker(P(d), 1/c*J(c)), P(c2))
    text <- paste(as.character(factor1), ":", as.character(factor3) )
  }
  if(H==6){
    K <- kronecker(kronecker(1/d*J(d), P(c)), P(c2))
    text <- paste( as.character(factor2), ":", as.character(factor3) )
  }
  # interaction effect of three factors
  if(H==7){
    K <- kronecker(kronecker(P(d), P(c)), P(c2))
    text <- paste(as.character(factor1), ":", as.character(factor2), ":", as.character(factor3) )
  }
  
  S <- 1
  # creating dual empirical covariance matrices
  K_B <- kronecker(S, K)
  V <- list(DualEmpirical2(Data = X, B=K)) #lapply(X, DualEmpirical2, B=K)
  
  ##########################
  ### U statistics
  #########################
  Q = data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU_onegroup(X,n,K)
    }
  }
  
  
  #################################################################################################
  
  # f
  f_1 <- 0
  f_2 <- 0
  
  for(i in 1:a){
    f_1 <- f_1 + (1*1/n[i])^2*.E1(n,i,V[[i]], nonparametric, Q)
    j <- i+1
    while(j<=a){
      f_1 <- f_1 + 2*(1*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j <- j+1
    }
  }
  
  for(i in 1:a){
    f_2 <- f_2 + (1*1/n[i])^2*.E2(n,i,V[[i]], nonparametric, Q)
    j <- i+1
    while(j<=a){
      f_2 <- f_2 + 2*S[i,j]*S[j,i]*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X,1/(n[j]-1)*K%*%t(X)%*%P(n[j])%*%X%*%K%*%t(X)%*%P(n[i]))
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
    f0_2 <- f0_2 + (1*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]], nonparametric, Q)
  }
  
  f0 <- f0_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit <- qf(1-alpha,f,f0)
  
  # Test
  
  direct <- 1/n[1]*var(X)
  eval.parent(substitute(varQGlobal <- direct))
  
  test <- (t(X_bar)%*%K_B%*%X_bar)/(t(rep(1,dim(K_B)[1]))%*%(K_B*direct)%*%(rep(1,dim(K_B)[1])))
  p.value <- 1-pf(test,f,f0)
  output <- data.frame(hypothesis=text,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# End ------------------------------------------------------------

