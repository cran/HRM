####################################################################################################################################
### Filename:    f2.R
### Description: Function for calculating the test statistic for one whole- and one subplot factor
###
###
###
####################################################################################################################################




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
hrm.1w.1f <- function(X, alpha, group , factor1, subject, data, H, text, nonparametric, ranked, varQGlobal ){

  stopifnot(is.data.frame(X),is.character(subject), is.character(group),is.character(factor1), alpha<=1, alpha>=0, is.logical(nonparametric))
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0


  group <- as.character(group)
  factor1 <- as.character(factor1)
  subject <- as.character(subject)

  X <- as.data.table(X)
  setnames(X, c(data, group, factor1, subject), c("data", "group", "factor1", "subject"))

  a <- nlevels(X[,group])
  d <- nlevels(X[,factor1])
  c <- 1
  n <- table(X[,group])/d

  if(nonparametric & is.null(ranked)) {
    X[,data:= 1/(sum(n)*d)*(pseudorank(X[,data], X[, group]) - 1/2)]
  }

  X <- split(X, X[,group], drop=TRUE)
  for(i in 1:a){
    X[[i]] <- X[[i]][ order(subject, factor1), ]
    X[[i]] <- X[[i]][,data]
    X[[i]] <- matrix(X[[i]],ncol=d*c,byrow=TRUE)
    n[i] <- dim(X[[i]])[1]
  }


  if(is.null(ranked)){
    eval.parent(substitute(ranked<-X))
  } else {
    X <- data.table::copy(ranked)
  }

  # creating X_bar (list with a entries)
  X_bar <- as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))

  if(H=="A"){
    K <- 1/d*J(d)
    S <- diag(n)-1/sum(n)*tcrossprod(n,n)
  } else if(H=="Au"){
    K <- 1/d*J(d)
    S <- P(a)
  } else if(H=="B"){
    K <- P(d)
    S <- J(a)
  } else if(H=="AB"){
    K <- P(d)
    S <- P(a)
  } else if(H=="A|B"){
    K <- I(d)
    S <- P(a)
  }

  # creating dual empirical covariance matrices
  K_AB <- kronecker(S, K)
  V <- lapply(X, DualEmpirical2, B=K)

  ##########################
  ### U statistics
  #########################

  Q = data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,K)
    }
  }



  #################################################################################################


  # f
  f_1 <- 0
  f_2 <- 0

  for(i in 1:a){
    f_1 <- f_1 + (S[i,i]*1/n[i])^2*.E1(n,i,V[[i]],nonparametric,Q)
    j <- i+1
    while(j<=a){
      f_1 <- f_1 + 2*(S[i,i]*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j <- j+1
    }
  }

  for(i in 1:a){
    f_2 <- f_2 + (S[i,i]*1/n[i])^2*.E2(n,i,V[[i]],nonparametric,Q)
    j <- i+1
    while(j<=a){
      f_2 <- f_2 + 2*S[i,j]*S[j,i]*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]], 1/(n[j]-1)*K%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%K%*%t(X[[i]])%*%P(n[i]))
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
    f0_2 <- f0_2 + (S[i,i]*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]],nonparametric,Q)
  }

  f0 <- f0_1/f0_2

  ##################################################################################################

  f <- abs(f)

  # critical value
  crit <- qf(1-alpha,f,f0)

  # variance estimator
  direct <- direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct <- direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }

  eval.parent(substitute(varQGlobal <- direct))

  den_one <- rep(1, dim(K_AB)[1])
  test <- crossprod(X_bar, crossprod(K_AB, X_bar))/crossprod(den_one, crossprod(K_AB*direct, den_one))
  p.value <- 1-pf(test,f,f0)
  output <- data.frame(hypothesis=text,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))

  return (output)
}

# End ------------------------------------------------------------
