####################################################################################################################################
### Filename:    f1.R
### Description: Function for calculating the test statistic for one subplot factor
###              
###
###
####################################################################################################################################

#' Test for no main effects and interactino effects of one between-subject factor and one crossed within-subject factors
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.1.one <- function(X, alpha , factor1, subject, data, formula ){

  temp0 <- hrm.1f(X, alpha , factor1,  subject, data, "B", paste(as.character(factor1)))

  output <- list()
  output$result <- rbind(temp0)
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c("none"), c(factor1))
  class(output) <- "HRM"
  
  return (output)
}


#' Test for interaction of factor A and B
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the response variable
#' @param H string specifying the hypothesis
#' @param text a string, which will be printed in the output
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.1f <- function(X, alpha , factor1, subject, data, H = "B", text ="" ){
  
  stopifnot(is.data.frame(X),is.character(subject), is.character(factor1), alpha<=1, alpha>=0)
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0  

  factor1 <- as.character(factor1)
  subject <- as.character(subject)
  a <- 1
  d <- nlevels(X[,factor1])
  c <- 1
  n <- rep(0,a) 
  
  for(i in 1:a){
    X<- X[ order(X[,subject], X[,factor1]), ]
    X<-X[,data]
    X<- matrix(X,ncol=d*c,byrow=TRUE)
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
  K_AB <- kronecker(S, K)
  V <- list(DualEmpirical2(Data = X, B=K)) #lapply(X, DualEmpirical2, B=K)
  
  #################################################################################################
  
  # f
  f_1 <- 0
  f_2 <- 0
  
  for(i in 1:a){
    f_1 <- f_1 + (1*1/n[i])^2*.E1(n,i,V[[i]])
    j <- i+1
    while(j<=a){
      f_1 <- f_1 + 2*(1*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j <- j+1
    }
  }
  
  for(i in 1:a){
    f_2 <- f_2 + (1*1/n[i])^2*.E2(n,i,V[[i]])
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
    f0_2 <- f0_2 + (1*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0 <- f0_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit <- qf(1-alpha,f,f0)
  
  # Test

  direct <- 1/n[1]*var(X)
  test <- (t(X_bar)%*%K_AB%*%X_bar)/(t(rep(1,dim(K_AB)[1]))%*%(K_AB*direct)%*%(rep(1,dim(K_AB)[1])))
  p.value <- 1-pf(test,f,f0)
  output <- data.frame(hypothesis=text,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# Hypothesis AC End ------------------------------------------------------------
