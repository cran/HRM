# functions which define the estimator

#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E1 = function(n,i, M) {
  return ((n[i]*(n[i]-1))/((n[i]-2)*(n[i]+1))*(matrix.trace(M)^2-2/n[i]*matrix.trace(M%*%M)))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E2 = function(n,i, M) {
  return ((n[i]-1)^2/((n[i]-2)*(n[i]+1))*(matrix.trace(M%*%M)-1/(n[i]-1)*matrix.trace(M)^2))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E3 = function(M_i, M_j) {
  return (matrix.trace(M_i)*matrix.trace(M_j))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E4 = function(M_i,M_j) {
  return (matrix.trace(M_i%*%M_j))
}

#######################################

#' Function for the output: significant p-values have on or more stars
#' 
#' @param value p-value
#' @export
#' @keywords internal
.hrm.sigcode = function(value) {
  
  char=""
  if(value <= 0.1 & value > 0.05) { char = "."}
  if(value <= 0.05 & value > 0.01) {char = '*'}
  if(value <= 0.01 & value > 0.001) {char = "**"}
  if(value <= 0.001 & value >= 0) {char = "***"}
  return (char)
}




#' Test for no main treatment effect (weighted version)
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_HAw.R
#' @export
hrm.A.weighted = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  
  f=0
  f0=0
  crit=0
  test=0
  
  # creating X_bar (list with a entries)
  X_bar = c(colMeans(X[[1]]))
  for(i in 2:a){
    X_bar = c(X_bar,colMeans(X[[i]]))
  }
  X_bar = as.matrix(X_bar)
  
  J_d = (rep(1,d)%*%t(rep(1,d)))
  D_a = diag(n)-1/sum(n)*n%*%t(n)
  K_A = kronecker(D_a,J_d)
  
  V = list(J_d%*%var(X[[1]])%*%t(J_d),J_d%*%var(X[[2]])%*%t(J_d))
  if(a>2){
    for(i in 3:a){
      V[[i]] = (J_d%*%var(X[[i]])%*%t(J_d))
    }
  }
  
  #################################################################################################
  # f = f_1/f_2: numerator degrees of freedom
  f_1 = 0
  f_2 = 0
  
  
  # computation of f_1
  for(i in 1:a){
    f_1 = f_1 + ((1-n[i]/sum(n))/d)^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*((1-n[i]/sum(n))/d)*((1-n[j]/sum(n))/d)*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  # computation of f_2
  for(i in 1:a){
    f_2 = f_2 + ((1-n[i]/sum(n))/d)^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2*(n[i]*n[j])/(d^2*sum(n)^2)*.E4(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  f=f_1/f_2
  
  
  ##################################################################################################
  
  
  
  #################################################################################################
  # f0 = f0_1/f0_2: denumerator degrees of freedom
  f0_1 = f_1
  f0_2 = 0
  
  # computation of f0_2
  for(i in 1:a){
    f0_2 = f0_2 + ((1-n[i]/sum(n))/d)^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0=f0_1/f0_2
  
  ##################################################################################################
  
  # critical value
  crit = qf(1-alpha,f,f0)
  
  # constructing the  Test
  direct = direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct = direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }
  test = (t(X_bar)%*%K_A%*%X_bar)/matrix.trace(K_A%*%direct)
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="A weighted",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  return (output)
}

# Hypothesis 1 ------------------------------------------------------------





#' Test for no main time effect
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_HB.R
#' @export
hrm.B = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  
  f=0
  f0=0
  crit=0
  test=0
  
  # creating X_bar (list with a entries)
  X_bar = c(colMeans(X[[1]]))
  for(i in 2:a){
    X_bar = c(X_bar,colMeans(X[[i]]))
  }
  X_bar = as.matrix(X_bar)
  
  J_a = (rep(1,a)%*%t(rep(1,a)))
  P_d = diag(d) - 1/d*(rep(1,d)%*%t(rep(1,d)))
  K_B = kronecker(J_a,P_d)
  
  V = list(P_d%*%var(X[[1]]),P_d%*%var(X[[2]]))
  if(a>2){
    for(i in 3:a){
      V[[i]] = (P_d%*%var(X[[i]]))
    }
  }
  
  #################################################################################################
  
  # f
  f_1 = 0
  f_2 = 0
  
  for(i in 1:a){
    f_1 = f_1 + (1/n[i])^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*(1/n[i])*(1/n[j])*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  for(i in 1:a){
    f_2 = f_2 + (1/n[i])^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2*1/(n[i]*n[j])*.E4(V[[i]],V[[j]])
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
    f0_2 = f0_2 + (1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]])
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
  test = (t(X_bar)%*%K_B%*%X_bar)/matrix.trace(K_B%*%direct)
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="B",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# Hypothesis 2 ------------------------------------------------------------







#' Test for no interaction between treatment and time
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_HAB.R
#' @export
hrm.AB = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  
  f=0
  f0=0
  crit=0
  test=0
  
  # creating X_bar (list with a entries)
  X_bar = c(colMeans(X[[1]]))
  for(i in 2:a){
    X_bar = c(X_bar,colMeans(X[[i]]))
  }
  X_bar = as.matrix(X_bar)
  
  P_a = diag(a) - 1/a*(rep(1,a)%*%t(rep(1,a)))
  P_d = diag(d) - 1/d*(rep(1,d)%*%t(rep(1,d)))
  K_AB = kronecker(P_a,P_d)
  
  V = list(P_d%*%var(X[[1]]),P_d%*%var(X[[2]]))
  if(a>2){
    for(i in 3:a){
      V[[i]] = (P_d%*%var(X[[i]]))
    }
  }
  
  #################################################################################################
  # f
  f_1 = 0
  f_2 = 0
  
  
  # phi = A
  for(i in 1:a){
    f_1 = f_1 + (1/n[i]*(1-1/a))^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*(1/n[i]*(1-1/a))*(1/n[j]*(1-1/a))*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  for(i in 1:a){
    f_2 = f_2 + (1/n[i]*(1-1/a))^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2/(n[i]*n[j]*a^2)*.E4(V[[i]],V[[j]])
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
    f0_2 = f0_2 + (1/n[i]*(1-1/a))^2*1/(n[i]-1)*.E2(n,i,V[[i]])
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
  test = (t(X_bar)%*%K_AB%*%X_bar)/matrix.trace(K_AB%*%direct)
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="AB",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# Hypothesis 3 ------------------------------------------------------------






#' Test for no simple treatment effect
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_HA_B.R
#' @export
hrm.A_B = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  
  f=0
  f0=0
  crit=0
  test=0
  
  # creating X_bar (list with a entries)
  X_bar = c(colMeans(X[[1]]))
  for(i in 2:a){
    X_bar = c(X_bar,colMeans(X[[i]]))
  }
  X_bar = as.matrix(X_bar)
  
  I_d = diag(d)
  P_a = diag(a) - 1/a*(rep(1,a)%*%t(rep(1,a)))
  K_A_B = kronecker(P_a,I_d)
  
  V = list(var(X[[1]]),var(X[[2]]))
  if(a>2){
    for(i in 3:a){
      V[[i]] = (var(X[[i]]))
    }
  }
  
  #################################################################################################
  # f
  f_1 = 0
  f_2 = 0
  
  for(i in 1:a){
    f_1 = f_1 + (1/n[i]*(1-1/a))^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*(1/n[i]*(1-1/a))*(1/n[j]*(1-1/a))*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  for(i in 1:a){
    f_2 = f_2 + (1/n[i]*(1-1/a))^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2/(n[i]*n[j]*a^2)*.E4(V[[i]],V[[j]])
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
    f0_2 = f0_2 + (1/n[i]*(1-1/a))^2*1/(n[i]-1)*.E2(n,i,V[[i]])
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
  test = (t(X_bar)%*%K_A_B%*%X_bar)/matrix.trace(K_A_B%*%direct)
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="A|B",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  return (output)
}

# Hypotheses 4 ------------------------------------------------------------




#' Test for no main treatment effect (unweighted version)
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_HAu.R
#' @export
hrm.A.unweighted = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  
  f=0
  f0=0
  crit=0
  test=0
  
  # creating X_bar (list with a entries)
  X_bar = c(colMeans(X[[1]]))
  for(i in 2:a){
    X_bar = c(X_bar,colMeans(X[[i]]))
  }
  X_bar = as.matrix(X_bar)
  
  P_a = diag(a) - 1/a*(rep(1,a)%*%t(rep(1,a)))
  J_d = (rep(1,d)%*%t(rep(1,d)))
  K_A = kronecker(P_a,J_d)
  
  V = list(J_d%*%var(X[[1]])%*%t(J_d),J_d%*%var(X[[2]])%*%t(J_d))
  if(a>2){
    for(i in 3:a){
      V[[i]] = (J_d%*%var(X[[i]])%*%t(J_d))
    }
  }
  
  #################################################################################################
  # f
  f_1 = 0
  f_2 = 0
  
  
  
  for(i in 1:a){
    f_1 = f_1 + ((1-1/a)/(d*n[i]))^2*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*((1-1/a)/(d*n[i]))*((1-1/a)/(d*n[j]))*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }
  
  for(i in 1:a){
    f_2 = f_2 + ((1-1/a)/(d*n[i]))^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2*(1/(a^2*n[i]*n[j]*d^2))*.E4(V[[i]],V[[j]])
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
    f0_2 = f0_2 + ((1-1/a)/(d*n[i]))^2*1/(n[i]-1)*.E2(n,i,V[[i]])
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
  test = (t(X_bar)%*%K_A%*%X_bar)/matrix.trace(K_A%*%direct)
  p.value=1-pf(test,f,f0)
  output = data.frame(hypothesis="A unweighted",df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(1-pf(test,f,f0)))
  

  
  return (output)
}

# Hypothesis 1/2 ------------------------------------------------------------


#' Test for no main treatment effect, no main time effect, no simple treatment effect and no interaction between treatment and time
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_all.R
#' @export
hrm.test = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  
  temp0=hrm.A.weighted(n,a,d,X,alpha)
  temp1 = hrm.A.unweighted(n,a,d,X,alpha)
  temp2=hrm.B(n,a,d,X,alpha)
  temp3=hrm.AB(n,a,d,X,alpha)
  temp4 = hrm.A_B(n,a,d,X,alpha)
  output = rbind(temp0,temp1,temp2,temp3,temp4)
  return (output)
}
