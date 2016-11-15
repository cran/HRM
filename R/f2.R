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
hrm.1w.1f = function(X, alpha, group , factor1, subject, data, H, text ){
  #Pa - Pd - 1/c Jc
  stopifnot(is.data.frame(X),is.character(subject), is.character(group),is.character(factor1), alpha<=1, alpha>=0)
  #stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  
  
  group = as.character(group)
  factor1 = as.character(factor1)
  subject = as.character(subject)
  X = split(X, X[,group], drop=TRUE)
  a = length(X)
  d = nlevels(X[[1]][,factor1])
  c = 1
  n = rep(0,a) 
  
  
  for(i in 1:a){
    X[[i]] = X[[i]][ order(X[[i]][,subject], X[[i]][,factor1]), ]
    X[[i]]=X[[i]][,data]
    X[[i]] = matrix(X[[i]],ncol=d*c,byrow=TRUE)
    n[i] = dim(X[[i]])[1]
  }
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  if(H=="A"){
    K = 1/d*J(d)
    S = diag(n)-1/sum(n)*n%*%t(n)
  }
  else if(H=="Au"){
    K = 1/d*J(d)
    S = P(a)
  }
  else if(H=="B"){
    K = P(d)
    S = J(a)
  }
  else if(H=="AB"){
    K = P(d)
    S = P(a)
  }
  
  # creating dual empirical covariance matrices
  K_AB = kronecker(S, K)
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
  test = (t(X_bar)%*%K_AB%*%X_bar)/(t(rep(1,dim(K_AB)[1]))%*%(K_AB*direct)%*%(rep(1,dim(K_AB)[1])))
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis=text,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# Hypothesis AC End ------------------------------------------------------------
