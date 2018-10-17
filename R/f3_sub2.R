####################################################################################################################################
### Filename:    f3_sub2.R
### Description: Function for calculating the test statistic for one whole- and two subplot factors
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
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param H string specifying the hypothesis
#' @param text a string, which will be printed in the output
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.1w.2f <- function(X, alpha, group , factor1, factor2, subject, data, H, text = "" , nonparametric, ranked, varQGlobal){

  stopifnot(is.data.frame(X),is.character(subject), is.character(group),is.character(factor1),is.character(factor2), alpha<=1, alpha>=0, is.logical(nonparametric))

  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0


  group <- as.character(group)
  factor1 <- as.character(factor1)
  factor2 <- as.character(factor2)
  subject <- as.character(subject)


  X <- as.data.table(X)
  setnames(X, c(data, group, factor1, factor2, subject), c("data", "group", "factor1", "factor2", "subject"))

  a <- nlevels(X[,group])
  d <- nlevels(X[,factor1])
  c <- nlevels(X[,factor2])
  n <- table(X[,group])/(d*c)

  if(nonparametric & is.null(ranked)) {
    X[,data:= 1/(sum(n)*d*c)*(pseudorank(X[,data], X[, group]) - 1/2)]
  }


  X <- split(X, X[,group], drop=TRUE)
  for(i in 1:a){
    X[[i]] <- X[[i]][ order(subject, factor1, factor2), ]
    X[[i]] <- X[[i]][,data]
    X[[i]] <- matrix(X[[i]],ncol=d*c,byrow=TRUE)
    n[i] <- dim(X[[i]])[1]
  }


  if(is.null(ranked)){
    eval.parent(substitute(ranked<-X))
  } else {
    X <- ranked
  }

  # creating X_bar (list with a entries)
  X_bar <- as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))

  # defining the hypothesis matrices
  if(H==1){ # A
    K <- 1/(d*c)*J(d*c)
    S <- P(a)
    text <- paste(as.character(group))
  } else if(H==2){ # B
    K <- kronecker(P(d), 1/c*J(c))
    S <- 1/a*J(a)
    text <- paste(as.character(factor1))
  } else if(H==3){ # C
    K <- kronecker(1/d*J(d), P(c))
    S <- 1/a*J(a)
    text <- paste(as.character(factor2))
  } else if(H==4){ # AB
    K <- kronecker(P(d), 1/c*J(c))
    S <- P(a)
    text <- paste(as.character(group),":",as.character(factor1))
  } else if(H==5){ # AC
    K <- kronecker(1/d*J(d), P(c))
    S <- P(a)
    text <- paste(as.character(group),":",as.character(factor2))
  } else if(H==6){ # BC
    K <- kronecker(P(d), P(c))
    S <- 1/a*J(a)
    text <- paste(as.character(factor1),":",as.character(factor2))
  } else if(H==7){ # ABC
    K <- kronecker(P(d), P(c))
    S <- P(a)
    text <- paste(as.character(group),":",as.character(factor1), ":", as.character(factor2))
  }

  # creating dual empirical covariance matrices
  K_Hypothesis <- kronecker(S, K)
  V <- lapply(X, DualEmpirical2, B=K)

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
    f_1 <- f_1 + (S[i,i]*1/n[i])^2*.E1(n,i,V[[i]], nonparametric, Q)
    j <- i+1
    while(j<=a){
      f_1 <- f_1 + 2*(S[i,i]*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j <- j+1
    }
  }

  for(i in 1:a){
    f_2 <- f_2 + (S[i,i]*1/n[i])^2*.E2(n,i,V[[i]], nonparametric, Q)
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
    f0_2 <- f0_2 + (S[i,i]*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]], nonparametric, Q)
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
  eval.parent(substitute(varQGlobal <- direct))

  test <- (t(X_bar)%*%K_Hypothesis%*%X_bar)/(t(rep(1,dim(K_Hypothesis)[1]))%*%(K_Hypothesis*direct)%*%(rep(1,dim(K_Hypothesis)[1])))
  p.value <- 1-pf(test,f,f0)
  output <- data.frame(hypothesis=text, df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))

  return (output)
}

# Hypothesis AC End ------------------------------------------------------------
