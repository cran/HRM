#' Test for 1 wholeplot and 3 subplot-factors
#' 
#' @param X dataframe containing the data and factor variables
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable crossed with the first one
#' @param factor3 column name of the data frame X of the third factor variable crossed with the first one
#' @param data column name of the data frame X containing the observed values
#' @param S Matrix for the wholeplot factor
#' @param K1 Matrix for the first subplot factor
#' @param K2 Matrix for the second subplot factor
#' @param K3 Matrix for the third subplot factor
#' @param hypothesis String which is printed in the console as ouput to indicate which hypothesis is tested.
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.1w.3f = function(X, alpha, group , factor1, factor2, factor3, subject, data, S, K1, K2, K3, hypothesis ){
  stopifnot(is.data.frame(X),is.character(subject), is.character(group),is.character(factor1),is.character(factor2), is.character(factor3),alpha<=1, alpha>=0)
  f=0
  f0=0
  crit=0
  test=0  
  
  # J_a - P_d - P_c
  
  group = as.character(group)
  factor1 = as.character(factor1)
  factor2 = as.character(factor2)
  subject = as.character(subject)
  X = split(X, X[,group], drop=TRUE)
  a = length(X)
  d = nlevels(X[[1]][,factor1])
  c = nlevels(X[[1]][,factor2])
  c2 =  nlevels(X[[1]][,factor3])
  n = rep(0,a) 
  
  
  for(i in 1:a){
    X[[i]] = X[[i]][ order(X[[i]][,subject], X[[i]][,factor1], X[[i]][,factor2]), ]
    X[[i]]=X[[i]][,data]
    X[[i]] = matrix(X[[i]],ncol=d*c*c2,byrow=TRUE)
    n[i] = dim(X[[i]])[1]
  }
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  
  # creating dual empirical covariance matrices
  S = S(a)
  K = kronecker(kronecker(K1(d), K2(c)), K3(c2))
  K_phi = kronecker(S, K)
  V = lapply(X, DualEmpirical2, B=K)
  
  #################################################################################################
  
  # f
  f_1 = 0
  f_2 = 0
  
  for(i in 1:a){
    f_1 = f_1 + (S[i,i]*1/n[i])^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*(S[i,i]*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  for(i in 1:a){
    f_2 = f_2 + (S[i,i]*1/n[i])^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2*S[i,j]*S[j,i]*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*K%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%K%*%t(X[[i]])%*%P(n[i]))
      j=j+1
    }
  }
  
  f=f_1/f_2
  
  
  ##################################################################################################
  
  
  
  #################################################################################################
  # f0
  f0_1 = f_1
  f0_2 = 0
  
  
  for(i in 1:a){
    f0_2 = f0_2 + (S[i,i]*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0=f0_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit = qf(1-alpha,f,f0)
  
  # Test
  
  direct = direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct = direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  test = (t(X_bar)%*%K_phi%*%X_bar)/(t(rep(1,dim(K_phi)[1]))%*%(K_phi*direct)%*%(rep(1,dim(K_phi)[1])))
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis=hypothesis,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# Hypothesis BC End ------------------------------------------------------------
