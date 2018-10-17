####################################################################################################################################
### Filename:    f1.R
### Description: Functions for calculating the test statistic for only one or two subplot factor
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
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.two <- function(X, alpha , factor1, factor2, subject, data, formula, testing = rep(1,3), nonparametric ){
  
  ranked <- NULL
  varQGlobal <- NULL
  
  temp0 <- if(testing[1]) { hrm.0w.2s(X, alpha , factor1,  factor2, subject, data, "B", paste(as.character(factor1)),nonparametric, ranked, varQGlobal) }
  temp1 <- if(testing[2]) { hrm.0w.2s(X, alpha , factor1,  factor2, subject, data, "C", paste(as.character(factor2)),nonparametric, ranked, varQGlobal) }
  temp2 <- if(testing[3]) { hrm.0w.2s(X, alpha , factor1,  factor2, subject, data, "BC", paste(as.character(factor1), ":", as.character(factor2) ),nonparametric, ranked, varQGlobal) }
  
  output <- list()
  output$result <- rbind(temp0,temp1,temp2)
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c("none"), c(factor1, factor2))
  output$data <- X
  output$var <- varQGlobal
  output$nonparametric <- nonparametric
  class(output) <- "HRM"
  
  return (output)
}



#' Test for one subplot factor
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.1.one <- function(X, alpha , factor1, subject, data, formula, nonparametric ){
  
  ranked <- NULL
  varQGlobal <- NULL
  
  temp0 <- hrm.1f(X, alpha , factor1,  subject, data, "B", paste(as.character(factor1)), nonparametric, ranked, varQGlobal)

  output <- list()
  output$result <- rbind(temp0)
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c("none"), c(factor1))
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
hrm.1f <- function(X, alpha , factor1, subject, data, H = "B", text ="" , nonparametric, ranked, varQGlobal){
  
  stopifnot(is.data.frame(X),is.character(subject), is.character(factor1), alpha<=1, alpha>=0, is.logical(nonparametric))
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0  
  
  factor1 <- as.character(factor1)
  subject <- as.character(subject)
  
  X <- as.data.table(X)
  setnames(X, c(data, factor1, subject), c("data", "factor1", "subject"))
  
  a <- 1
  d <- nlevels(X[,factor1])
  c <- 1
  n <- dim(X)[1]

  if(nonparametric){
    X[,data := (rank(X[,data], ties.method = "average")-1/2)*1/(n[1]*a*d)]
  }
  
  for(i in 1:a){
    X <- X[ order(subject, factor1), ]
    X <- X[,data]
    X <- matrix(X,ncol=d*c,byrow=TRUE)
    n[i] <- dim(X)[1]
  }

  # creating X_bar (list with a entries)
  X_bar <- colMeans(X)  # as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  
  


 if(H=="B"){
    K <- P(d)
    text <- paste(as.character(factor1))
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
hrm.0w.2s <- function(X, alpha , factor1, factor2, subject, data, H = "B", text ="", nonparametric, ranked, varQGlobal ){
  
  stopifnot(is.data.frame(X),is.character(subject), is.character(factor1), is.character(factor2), alpha<=1, alpha>=0)
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0  
 
  factor1 <- as.character(factor1)
  factor2 <- as.character(factor2)
  subject <- as.character(subject)
  
  X <- as.data.table(X)
  setnames(X, c(data, factor1, subject, factor2), c("data", "factor1", "subject", "factor2"))
  
  a <- 1
  d <- nlevels(X[,factor1])
  c <- nlevels(X[,factor2])
  n <- dim(X)[1]
  
  if(nonparametric & is.null(ranked)){ 
    X[,data := (rank(X[,data], ties.method = "average")-1/2)*1/(n[1]*a*d*c)]
    #X[,data]<- (rank(X[,data], ties.method = "average")-1/2)*1/(n[1]*a*d*c)
  }

  for(i in 1:a){
    X <- X[ order(subject, factor1, factor2), ]
    X <- X[,data]
    X <- matrix(X,ncol=d*c,byrow=TRUE)
    n[i] <- dim(X)[1]
  }
  
  if(is.null(ranked)){
    eval.parent(substitute(ranked<-X))
  } else {
    X <- ranked
  }

  # creating X_bar (list with a entries)
  X_bar <- colMeans(X)  # as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  
  
  
  
  if(H=="B"){
    K <- kronecker(P(d), 1/c*J(c))
  }
  if(H=="C"){
    K <- kronecker(1/d*J(d), P(c))
  }
  if(H=="BC"){
    K <- kronecker(P(d), P(c))
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

