#' EEG data of 160 subjects
#' 
#' A dataset containing EEG data (Staffen et al., 2014) of 160 subjects, 4 variables are measured at ten different locations.
#' The columns are as follows:
#' 
#' \itemize{
#'   \item group. Diagnostic group of the subject: Alzheimer's Disease (AD), Mild Cognitive Impairment (MCI), Subject Cognitive Complaints (SCC+, SCC-).
#'   \item value. Measured data of a subject at a specific variable and region.
#'   \item sex. Sex of the subject: Male (M) or Female (W).
#'   \item subject. A unique identification of a subject.
#'   \item variable. The variales measured are activity, complexity, mobility and brain rate coded from 1 to 4.
#'   \item region. Frontal left/right, central left/right, temporal left/right, occipital left/right, parietal left/right coded as 1 to 10. 
#'   \item dimension. Mixing variable and region together, levels range from 1 to 40.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name EEG
#' @usage data(EEG)
#' @format A data frame with 6400 rows and 7 variables.
"EEG"

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
#' @keywords internal
.hrm.sigcode = function(value) {
  
  char=""
  if(value <= 0.1 & value > 0.05) { char = "."}
  if(value <= 0.05 & value > 0.01) {char = '*'}
  if(value <= 0.01 & value > 0.001) {char = "**"}
  if(value <= 0.001 & value >= 0) {char = "***"}
  return (char)
}

I = function(size){
  return (diag(rep(1,size)))
}

J = function(size){
  return (rep(1,size)%*%t(rep(1,size)))
}

P = function(size){
  return (I(size)-J(size)*1/size)
}


DualEmpirical = function(Data, B){
  n = dim(Data)[1]
  B = B(dim(Data)[2]) # B = e.g. t(J_d)%*%J_d
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}


DualEmpirical2 = function(Data, B){
  n = dim(Data)[1]
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}

#' Test for no main treatment effect (weighted version)
#' 
#' @param n an vector containing the sample sizes of all groups
#' @param a number of groups
#' @param d number of dimensions (time points)
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.A.weighted = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0
  if(is.data.frame(X[[1]])){
    X = lapply(X, as.matrix)
  }
  
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  D_a = diag(n)-1/sum(n)*n%*%t(n)
  K_A = kronecker(D_a,J(d))
  
  # creating dual empirical covariance matrices
  V = lapply(X, DualEmpirical, B = J)
  
  

  #################################################################################################
  # f = f_1/f_2: numerator degrees of freedom
  f_1 = 0
  f_2 = 0

  # computation of f_1
  for(i in 1:a){
    f_1 = f_1 + ((1-n[i]/sum(n))^2)*.E1(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_1 = f_1 + 2*((1-n[i]/sum(n)))*((1-n[j]/sum(n)))*.E3(V[[i]],V[[j]])
      j=j+1
    }
  }

  # computation of f_2
  for(i in 1:a){
    f_2 = f_2 + ((1-n[i]/sum(n)))^2*.E2(n,i,V[[i]])
    j=i+1
    while(j<=a){
      f_2 = f_2 + 2*(n[i]*n[j])/(d^2*sum(n)^2)*.E4(d/(n[i]-1)*P(n[i])%*%X[[i]],d/(n[j]-1)*J(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%J(d)%*%t(X[[i]])%*%P(n[i]))
      j=j+1
    }
  }
  
  f=f_1/f_2

  ##################################################################################################
  
  
  
  #################################################################################################
  # f0 = f0_1/f0_2: denumerator degrees of freedom
  # f0_1 = f_1 numerator are the same
  f0_2 = 0
  
  # computation of f0_2
  for(i in 1:a){
    f0_2 = f0_2 + ((1-n[i]/sum(n)))^2*1/(n[i]-1)*.E2(n,i,V[[i]])
  }
  
  f0=f_1/f0_2
  
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

  test = (t(X_bar)%*%K_A%*%X_bar)/(t(rep(1,dim(K_A)[1]))%*%(K_A*direct)%*%(rep(1,dim(K_A)[1])))
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
#' @keywords internal
hrm.B = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  X = lapply(X, as.matrix)
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  D_a = diag(n)-1/sum(n)*n%*%t(n)
  K_A = kronecker(D_a,J(d))
  
  # creating dual empirical covariance matrices
  V = lapply(X, DualEmpirical, B = P)
  
  K_B = kronecker(J(a),P(d))
  
  #   V = list(P_d%*%var(X[[1]]),P_d%*%var(X[[2]]))
  #   if(a>2){
  #     for(i in 3:a){
  #       V[[i]] = (P_d%*%var(X[[i]]))
  #     }
  #   }
  
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
      f_2 = f_2 + 2*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*P(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%P(d)%*%t(X[[i]])%*%P(n[i]))
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
  test = (t(X_bar)%*%K_B%*%X_bar)/(t(rep(1,dim(K_B)[1]))%*%(K_B*direct)%*%(rep(1,dim(K_B)[1])))
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
#' @keywords internal
hrm.AB = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  X = lapply(X, as.matrix)
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  # creating dual empirical covariance matrices
  V = lapply(X, DualEmpirical, B = P)
  
  K_AB = kronecker(P(a),P(d))
  
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
      f_2 = f_2 + 2/(n[i]*n[j]*a^2)*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*P(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%P(d)%*%t(X[[i]])%*%P(n[i]))
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
  test = (t(X_bar)%*%K_AB%*%X_bar)/(t(rep(1,dim(K_AB)[1]))%*%(K_AB*direct)%*%(rep(1,dim(K_AB)[1])))
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
#' @keywords internal
hrm.A_B = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  X = lapply(X, as.matrix)
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  # creating dual empirical covariance matrices
  V = lapply(X, DualEmpirical, B = I)
  
  K_A_B = kronecker(P(a),I(d))
  
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
      f_2 = f_2 + 2/(n[i]*n[j]*a^2)*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*t(X[[j]])%*%P(n[j])%*%X[[j]]%*%t(X[[i]])%*%P(n[i]))
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
  test = (t(X_bar)%*%K_A_B%*%X_bar)/(t(rep(1,dim(K_A_B)[1]))%*%(K_A_B*direct)%*%(rep(1,dim(K_A_B)[1])))
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
#' @keywords internal
hrm.A.unweighted = function(n, a, d, X, alpha){
  
  stopifnot(is.list(X), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  stopifnot(a == length(X))
  f=0
  f0=0
  crit=0
  test=0  
  X = lapply(X, as.matrix)
  
  # creating X_bar (list with a entries)
  X_bar = as.matrix(vec(sapply(X, colMeans, na.rm=TRUE)))
  
  # creating dual empirical covariance matrices, where observations are transformed by a matrix B
  V = lapply(X, DualEmpirical, B = J)
  
  K_A = kronecker(P(a),J(d))
  

  
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
      f_2 = f_2 + 2*(1/(a^2*n[i]*n[j]*d^2))*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*J(d)%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%J(d)%*%t(X[[i]])%*%P(n[i]))
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
  test = (t(X_bar)%*%K_A%*%X_bar)/(t(rep(1,dim(K_A)[1]))%*%(K_A*direct)%*%(rep(1,dim(K_A)[1])))
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
#' @param data A list containing the data matrices of all groups. The rows are the independent subjects, these observations are assumed to be multivariate normally distributed. The columsn of all matrices need to be in the same order.
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_1.txt
#' @keywords internal
hrm.test = function(n, a, d, data, alpha=0.05){
  .Deprecated("hrm.test.matrix", package=NULL,
              old = as.character(sys.call(sys.parent()))[1L])
  
  return(hrm.test.matrix(data,alpha))
}

#' Test for no main treatment effect, no main time effect, no simple treatment effect and no interaction between treatment and time
#' 
#' @param data A list containing the data matrices of all groups. The rows are the independent subjects, these observations are assumed to be multivariate normally distributed. The columsn of all matrices need to be in the same order.
#' @param alpha alpha level used for the test
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_1.txt
#' @keywords export
hrm.test.matrix = function(data, alpha=0.05){
  
  if(!is.list(data)){
    stop("data needs to be a list containing the data matrices of all groups")
  }
  
  a = length(data)
  if(a < 2){
    stop("At least two groups are needed.")
  }
  
  n = rep(0,a)
  for(i in 1:a){
    tmp = data[[i]]
    if(!is.matrix(tmp)){
      stop("The elements of data need to be matrices.")
    }
    n[i] = dim(tmp)[1]
  }
  d = dim(data[[1]])[2]
  d2 = dim(data[[2]])[2]
  if(d != d2){
    stop("The number of measurements for each group need to be the same.")
  }
  
  stopifnot(is.list(data), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  
  temp0=hrm.A.weighted(n,a,d,data,alpha)
  temp1 = hrm.A.unweighted(n,a,d,data,alpha)
  temp2=hrm.B(n,a,d,data,alpha)
  temp3=hrm.AB(n,a,d,data,alpha)
  temp4 = hrm.A_B(n,a,d,data,alpha)
  output = rbind(temp0,temp1,temp2,temp3,temp4)
  return (output)
}

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
hrm.test.2.one = function(X, alpha, group , factor1, subject, data, testing = rep(1,4) ){
  
  temp0= if(testing[1]) {hrm.1w.1f(X, alpha, group , factor1,  subject, data, "A", paste(as.character(group), " (weighted)"))}
  temp1= if(testing[2]) {hrm.1w.1f(X, alpha, group , factor1,  subject, data, "Au", paste(as.character(group)))}
  temp2= if(testing[3]) {hrm.1w.1f(X, alpha, group , factor1,  subject, data, "B", paste(as.character(factor1)))}
  temp3= if(testing[4]) {hrm.1w.1f(X, alpha, group , factor1, subject, data, "AB", paste(as.character(group), ":",as.character(factor1)))}

  output = rbind(temp0, temp1, temp2, temp3)
  return (output)
}

#' Test for no main effects and interactino effects of one between-subject factor and two crossed within-subject factors
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.within = function(X, alpha, group , factor1, factor2, subject, data, testing = rep(1,7) ){
  
  temp0= if(testing[1]) {hrm.1w.2f.A(X, alpha, group , factor1, factor2, subject, data )}
  temp1= if(testing[2]) {hrm.1w.2f.B(X, alpha, group , factor1, factor2, subject, data )}
  temp2= if(testing[3]) {hrm.1w.2f.C(X, alpha, group , factor1, factor2, subject, data )}
  temp3= if(testing[4]) {hrm.1w.2f.AB(X, alpha, group , factor1, factor2, subject, data )}
  temp4= if(testing[5]) {hrm.1w.2f.AC(X, alpha, group , factor1, factor2, subject, data )}
  temp5= if(testing[6]) {hrm.1w.2f.BC(X, alpha, group , factor1, factor2, subject, data )}
  temp6= if(testing[7]) {hrm.1w.2f.ABC(X, alpha, group , factor1, factor2, subject, data )}
  
  output = rbind(temp0, temp1, temp2, temp3,temp4,temp5,temp6)
  return (output)
}



#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor column name of the data frame X of within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.between = function(X, alpha, group , subgroup, factor, subject, data, testing = rep(1,7) ){
  
  temp1 = if(testing[1]) {hrm.2w.1f.A(X, alpha, group , subgroup, factor, subject, data )}
  temp2 = if(testing[2]) {hrm.2w.1f.A2(X, alpha, group , subgroup, factor, subject, data )}
  temp3 = if(testing[3]) {hrm.2w.1f.B(X, alpha, group , subgroup, factor, subject, data  )}
  temp4 = if(testing[4]) {hrm.2w.1f.AA2(X, alpha, group , subgroup, factor, subject, data  )}
  temp5 = if(testing[5]) {hrm.2w.1f.AB(X, alpha, group , subgroup, factor, subject, data  )}
  temp6 = if(testing[6]) {hrm.2w.1f.A2B(X, alpha, group , subgroup, factor, subject, data  )}
  temp7 = if(testing[7]) {hrm.2w.1f.AA2B(X, alpha, group , subgroup, factor, subject, data  )}
  
  output = rbind(temp1, temp2, temp3, temp4, temp5, temp6, temp7)
  return (output)
}


#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor1 column name of the data frame X of the first within-subject factor
#' @param factor2 column name of the data frame X of the second within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.between.within = function(X, alpha, group , subgroup, factor1, factor2, subject, data, testing = rep(1,15) ){
  temp0= if(testing[1]) {hrm.2w.2f.A(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp1= if(testing[2]) {hrm.2w.2f.A2(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp2= if(testing[3]) {hrm.2w.2f.B(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp3= if(testing[4]) {hrm.2w.2f.C(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp4= if(testing[5]) {hrm.2w.2f.AA2(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp5= if(testing[6]) {hrm.2w.2f.AB(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp6= if(testing[7]) {hrm.2w.2f.AC(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp7= if(testing[8]) {hrm.2w.2f.A2B(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp8= if(testing[9]) {hrm.2w.2f.A2C(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp9= if(testing[10]) {hrm.2w.2f.BC(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp10= if(testing[11]) {hrm.2w.2f.AA2B(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp11= if(testing[12]) {hrm.2w.2f.AA2C(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp12= if(testing[13]) {hrm.2w.2f.ABC(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp13= if(testing[14]) {hrm.2w.2f.A2BC(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  temp14= if(testing[15]) {hrm.2w.2f.AA2BC(X, alpha, group , subgroup, factor1, factor2, subject, data )}
  
  
  output = rbind(temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14)
  return (output)
}


#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#' 
#' @param X list containing the data matrices of all groups
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first within-subject factor
#' @param factor2 column name of the data frame X of the second within-subject factor
#' @param factor3 column name of the data frame X of the third within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.3.between = function(X, alpha, group , factor1, factor2, factor3, subject, data, testing = rep(1,15) ){
  temp0= if(testing[1]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3, subject, data, P, J, J, J,  paste(as.character(group) ) )}
  temp1= if(testing[2]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, P, J, J, paste(as.character(factor1)) )}
  temp2= if(testing[3]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, J, P, J, paste(as.character(factor2)) )}
  temp3= if(testing[4]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, J, J, P, paste(as.character(factor3)) )}
  temp4= if(testing[5]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, P, J, J, paste(as.character(group),":",as.character(factor1)) )} 
  temp5= if(testing[6]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, J, P, J, paste(as.character(group),":",as.character(factor2)) )}
  temp6= if(testing[7]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, J, J, P, paste(as.character(group),":",as.character(factor3)) )}
  temp7= if(testing[8]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, P, P, J, paste(as.character(factor1),":",as.character(factor2)) )}
  temp8= if(testing[9]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, P, J, P, paste(as.character(factor1),":",as.character(factor3)) )}
  temp9= if(testing[10]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, J, P, P,paste(as.character(factor2),":",as.character(factor3)) )}
  temp10= if(testing[11]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, P, P, J, paste(as.character(group),":",as.character(factor1), ":", as.character(factor2)) )}
  temp11= if(testing[12]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, P, J, P, paste(as.character(group),":",as.character(factor1), ":", as.character(factor3)) )}
  temp12= if(testing[13]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, J, P, P, paste(as.character(group),":",as.character(factor2), ":", as.character(factor3)) )}
  temp13= if(testing[14]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, P, P, P, paste(as.character(factor1),":",as.character(factor2), ":", as.character(factor3)) )}
  temp14= if(testing[15]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, P, P, P, paste(as.character(group),":",as.character(factor1), ":", as.character(factor2), ":", as.character(factor3)) )}
  
  
  output = rbind(temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14)
  return (output)
}


#' Test for main effects and interaction effects of one or two between-subject factors and one or two within-subject factors
#' 
#' @param X A data.frame containing the data
#' @param alpha alpha level used for the test
#' @param group column name within the data frame data specifying the groups
#' @param subgroup column name within the data frame data specifying the subgroups (crossed with groups)
#' @param factor1 column name within the data frame data specifying the first subplot-factor
#' @param factor2 column name within the data frame data specifying the the second subplot-factor (crossed with factor1)
#' @param subject column name within the data frame X identifying the subjects
#' @param data column name within the data frame X containing the response variable
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_2.txt
#' @keywords internal
hrm.test.2 = function(X, alpha = 0.05, group , subgroup, factor1, factor2, subject, data ){
  .Deprecated("hrm.test.dataframe", package=NULL,
              old = as.character(sys.call(sys.parent()))[1L])
  return(hrm.test.dataframe(X, alpha, group , subgroup, factor1, factor2, subject, data ))
}

#' Test for main effects and interaction effects of one or two between-subject factors and one, two or three within-subject factors (at most four factors can be used).
#' 
#' @param data A data.frame containing the data
#' @param alpha alpha level used for the test
#' @param group column name within the data frame data specifying the groups
#' @param subgroup column name within the data frame data specifying the subgroups (crossed with groups)
#' @param factor1 column name within the data frame data specifying the first subplot-factor
#' @param factor2 column name within the data frame data specifying the the second subplot-factor (crossed with factor1)
#' @param subject column name within the data frame X identifying the subjects
#' @param response column name within the data frame X containing the response variable
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_2.txt
#' @keywords export
hrm.test.dataframe = function(data, alpha = 0.05, group , subgroup, factor1, factor2, subject, response ){
  
  temp = data
  data = response
  X = temp
  
  if(missing(group) || !is.character(group)){
    print("At least one between-subject factor is needed!")
    stop("group column name not specified ")
  }
  
  if(missing(factor1) || !is.character(factor1)){
    print("At least one between-subject factor is needed!")
    stop("factor1 column name not specified ")
  }
  
  if(missing(X) || !is.data.frame(X)){
    stop("dataframe needed")
  }
  
  if(missing(subject) || !is.character(subject)){
    stop("subject column name not specified")
  }
  
  if(missing(data) || !is.character(data)){
    stop("data column name not specified")
  }
  
  if(sum(is.na(X[,data]))>=1){
    warning("Your data contains missing values!")
  }
  
  if(!is.double(alpha)){
    stop("alpha level needs to be a number between 0 and 1")
  }
  if(is.double(alpha)){
    if(alpha > 1 || alpha < 0){
      stop("alpha level needs to be a number between 0 and 1")
    }
  }
  dat = X
  dat = data.frame(dat, subj = X[,subject])
  s1=subset(dat, dat$subj==dat$subj[1])
  measurements = dim(s1)[1]
  countSubplotFactor = nlevels(s1[,factor1])
  if(!missing(factor2)){
    countSubplotFactor = countSubplotFactor*nlevels(s1[,factor2])
  }
  if(!(measurements == countSubplotFactor)){
    stop(paste("The number of repeated measurements per subject (", measurements, ") is uneqal to the number of levels of the subplot factors (", countSubplotFactor, ")."))
  }
  if(missing(factor2) & missing(subgroup)){
    if(!is.factor(X[,group])){
      stop(paste("The column ", group, " is not a factor." ))
    }
    if(!is.factor(X[,factor1])){
      stop(paste("The column ", factor1, " is not a factor." ))
    }
    return(hrm.test.2.one(X, alpha, group , factor1, subject, data ))
  }
  
  if(missing(factor2) & !missing(subgroup)){
    if(!is.character(subgroup)){
      stop("subgroup column name not specified")
    }
    if(!is.factor(X[,group])){
      stop(paste("The column ", group, " is not a factor." ))
    }
    if(!is.factor(X[,subgroup])){
      stop(paste("The column ", subgroup, " is not a factor." ))
    }
    if(!is.factor(X[,factor1])){
      stop(paste("The column ", factor1, " is not a factor." ))
    }
    return(hrm.test.2.between(X, alpha, group , subgroup, factor1, subject, data ))
  }
  
  if(!missing(factor2) & missing(subgroup)){
    if(!is.character(factor2)){
      stop("factor2 column name not specified")
    }
    if(!is.factor(X[,group])){
      stop(paste("The column ", group, " is not a factor." ))
    }
    if(!is.factor(X[,factor2])){
      stop(paste("The column ", factor2, " is not a factor." ))
    }
    if(!is.factor(X[,factor1])){
      stop(paste("The column ", factor1, " is not a factor." ))
    }
    return(hrm.test.2.within(X, alpha, group , factor1, factor2, subject, data ))
  }
  
  if(!missing(factor2) & !missing(subgroup)){
    if(!is.character(factor2)){
      stop("factor2 column name not specified")
    }
    if(!is.character(subgroup)){
      stop("subgroup column name not specified")
    }
    if(!is.factor(X[,group])){
      stop(paste("The column ", group, " is not a factor." ))
    }
    if(!is.factor(X[,subgroup])){
      stop(paste("The column ", subgroup, " is not a factor." ))
    }
    if(!is.factor(X[,factor1])){
      stop(paste("The column ", factor1, " is not a factor." ))
    }
    if(!is.factor(X[,factor2])){
      stop(paste("The column ", factor2, " is not a factor." ))
    }
    return(hrm.test.2.between.within(X, alpha, group , subgroup, factor1, factor2, subject, data ))
  }
  
}



#' Test for main effects and interaction effects of one or two between-subject factors and one, two or three within-subject factors (at most four factors can be used).
#' 
#' @param data A data.frame containing the data. The columns containing the factor variables need to have the type 'factor'. One column is needed to indentify the subjects.
#' @param alpha alpha level used for the test
#' @param formula A model formula object. The left hand side contains the response variable and the right hand side contains the whole- and subplot factors.
#' @param subject column name within the data frame X identifying the subjects
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_3.txt
#' @keywords export
hrm_test = function(formula, data, alpha = 0.05,  subject ){
  
  if(missing(data) || !is.data.frame(data)){
    stop("dataframe needed")
  }
  if(missing(subject) || !is.character(subject)){
    stop("subject column name not specified")
  }
  if(!is.double(alpha)){
    stop("alpha level needs to be a number between 0 and 1")
  }
  if(is.double(alpha)){
    if(alpha > 1 || alpha < 0){
      stop("alpha level needs to be a number between 0 and 1")
    }
  }
  dat = model.frame(formula, data)
  dat2 <- data.frame(dat,subj=data[,subject])
  m = ncol(dat)
  
  if(!is.numeric(dat[,1])){
    stop("Response variable needs to be numeric!")
  }
  
  
  
  # find out, in which columns are the wholeplot or subplot factors
  s1=subset(dat2, dat2$subj==dat2$subj[1])
  measurements = dim(s1)[1]
  countSubplotFactor = 1
  wholeplot=rep(-1, m)
  subplot=rep(-1, m)
  for(i in 2:m){
    if(!is.factor(dat2[,i])){
      stop(paste("The column ", colnames(dat2)[i], " is not a factor." ))
    }
    if(length(unique(s1[,i]))==nlevels(dat2[,i])){
      subplot[i]=1
      countSubplotFactor = countSubplotFactor*nlevels(s1[,i])
    }
    else{
      wholeplot[i]=1
    }
  }
  
  wholeplot = which(wholeplot==1)
  subplot = which( subplot==1)
  if(!(measurements == countSubplotFactor)){
    stop(paste("The number of repeated measurements per subject (", measurements, ") is uneqal to the number of levels of the subplot factors (", countSubplotFactor, ")."))
  }
  if(length(wholeplot)>2){
    stop("Too many factors are used! Only two wholelot-factors are supported.")
  }
  if(length(subplot)>3){
    stop("Too many factors are used! Only three subplotlot-factors are supported.")
  }
  if(length(wholeplot)>1 & length(subplot)==3){
    stop("Too many factors are used! Only one whole- and three subplot-factors are supported.")
  }
  if(length(wholeplot)<1 || length(subplot)<1){
    stop("The model needs at least one between-subject factors and one within-subject factors.")
  }
  
  if(length(wholeplot)==1 & length(subplot)==1){
    group = colnames(dat2)[wholeplot[1]]
    factor1 = colnames(dat2)[subplot[1]]
    x=attributes(terms.formula(formula))$term.labels
    
    wholeplot=colnames(dat2[,wholeplot])
    subplot=colnames(dat2[,subplot])
    
    testing = rep(0,4)
    for(i in 1:length(x)){
      
      tmp = strsplit(x[i],":")
      l = length(tmp[[1]])
      
      # interaction hypothesis of 2 factors
      if(l == 2){testing[4]=1}
      
      # l = 1
      else{
        if(group == x[i]){
          testing[1]=1
          testing[2]=1
        }
        else if(factor1 == x[i]){testing[3]=1}
      }
    }
    X=data
    data = colnames(dat)[1]
    return(hrm.test.2.one(X, alpha, group , factor1, subject, data, testing ))   
  }
  
  if(length(wholeplot)==2 & length(subplot)==1){
    group = colnames(dat2)[wholeplot[1]]
    subgroup = colnames(dat2)[wholeplot[2]]
    factor1 = colnames(dat2)[subplot[1]]
    x=attributes(terms.formula(formula))$term.labels
    
    wholeplot=colnames(dat2[,wholeplot])
    subplot=colnames(dat2[,subplot])
    
    testing = rep(0,7)
    for(i in 1:length(x)){
      
      tmp = strsplit(x[i],":")
      l = length(tmp[[1]])
      
      # interaction hypothesis of 4 factors
      if(l == 3){testing[7]=1}
      
      # find out which interaction hypothesis of 2 factors is tested
      else if(l==2){
        if(grepl(group,x[i])){
          if(grepl(subgroup,x[i])){
            testing[4]=1
          }
          if(grepl(factor1,x[i])){
            testing[5]=1
          }
        }
        else if(grepl(subgroup,x[i])){
          if(grepl(factor1,x[i])){
            testing[6]=1
          }
        }
      }
      # l = 1
      else{
        if(group == x[i]){testing[1]=1}
        else if(subgroup == x[i]){testing[2]=1}
        else if(factor1 == x[i]){testing[3]=1}
      }
    }
    X=data
    data = colnames(dat)[1]
    return(hrm.test.2.between(X, alpha, group , subgroup, factor1, subject, data, testing )) 
  }
  
  if(length(wholeplot)==1 & length(subplot)==2){
    group = colnames(dat2)[wholeplot[1]]
    factor1 = colnames(dat2)[subplot[1]]
    factor2 = colnames(dat2)[subplot[2]]
    x=attributes(terms.formula(formula))$term.labels
    
    wholeplot=colnames(dat2[,wholeplot])
    subplot=colnames(dat2[,subplot])
    
    testing = rep(0,7)
    for(i in 1:length(x)){
      
      tmp = strsplit(x[i],":")
      l = length(tmp[[1]])
      
      # interaction hypothesis of 3 factors
      if(l == 3){testing[7]=1}
      
      # find out which interaction hypothesis of 2 factors is tested
      else if(l==2){
        if(grepl(group,x[i])){
          if(grepl(factor1,x[i])){
            testing[4]=1
          }
          if(grepl(factor2,x[i])){
            testing[5]=1
          }
        }
        else if(grepl(factor1,x[i])){
          testing[6]=1
        }
      }
      # l = 1
      else{
        if(group == x[i]){testing[1]=1}
        else if(factor1 == x[i]){testing[2]=1}
        else {testing[3]=1}
      }
    }
    X=data
    data = colnames(dat)[1]
    return(hrm.test.2.within(X, alpha, group , factor1, factor2, subject, data, testing ))
  }
  
  if(length(wholeplot)==1 & length(subplot)==3){
    group = colnames(dat2)[wholeplot[1]]
    factor1 = colnames(dat2)[subplot[1]]
    factor2 = colnames(dat2)[subplot[2]]
    factor3 = colnames(dat2)[subplot[3]]
    
    x=attributes(terms.formula(formula))$term.labels
    
    wholeplot=colnames(dat2[,wholeplot])
    subplot=colnames(dat2[,subplot])
    
    testing = rep(0,15)
    for(i in 1:length(x)){
      
      tmp = strsplit(x[i],":")
      l = length(tmp[[1]])
      
      # interaction hypothesis of 3 factors
      if(l == 4){testing[15]=1}
      
      # find out which interaction hypothesis of 3 factors is tested
      else if(l==3){
        if(grepl(group,x[i])){
          if(grepl(factor1,x[i])){
            if(grepl(factor2,x[i])){
              testing[11]=1
            }
            else{
              testing[12]=1
            }
          }
          if(grepl(factor2,x[i])){
            if(grepl(factor3,x[i])){
              testing[13]=1
            }
          }
        }
        else if(grepl(factor1,x[i])){
          testing[14]=1
        }
      }
      # l = 2
      else if (l==2){
        if(grepl(group,x[i])){
          if(grepl(factor1,x[i])){
            testing[5]=1
          }
          if(grepl(factor2,x[i])){
            testing[6]=1
          }
          if(grepl(factor3,x[i])){
            testing[7]=1
          }
        }
        else if(grepl(factor1,x[i])){
          if(grepl(factor2,x[i])){
            testing[8]=1
          }
          if(grepl(factor3,x[i])){
            testing[9]=1
          }
        }
        else if(grepl(factor2,x[i])){
          if(grepl(factor3,x[i])){
            testing[10]=1
          }
        }
      }
      else if(l==1){
        if(grepl(group, x[[i]])){
          testing[1]=1
        }
        if(grepl(factor1, x[[i]])){
          testing[2]=1
        }
        if(grepl(factor2, x[[i]])){
          testing[3]=1
        }
        if(grepl(factor3, x[[i]])){
          testing[4]=1
        }
      }
    }
    X=data
    data = colnames(dat)[1]
    return(hrm.test.3.between(X, alpha, group , factor1, factor2, factor3, subject, data, testing ))
  }
  
  
  
  if(length(wholeplot)==2 & length(subplot)==2){
    group = colnames(dat2)[wholeplot[1]]
    subgroup = colnames(dat2)[wholeplot[2]]
    factor1 = colnames(dat2)[subplot[1]]
    factor2 = colnames(dat2)[subplot[2]]
    x=attributes(terms.formula(formula))$term.labels
    
    wholeplot=colnames(dat2[,wholeplot])
    subplot=colnames(dat2[,subplot])
    
    testing = rep(0,15)
    for(i in 1:length(x)){
      
      tmp = strsplit(x[i],":")
      l = length(tmp[[1]])
      
      # interaction hypothesis of 4 factors
      if(l == 4){testing[15]=1}
      
      # find out which interaction hypothesis of 3 factors is tested
      else if(l==3){
        if(grepl(group,x[i])){
          if(grepl(subgroup,x[i])){
            if(grepl(factor1,x[i])){
              testing[11]=1
            }
            else{
              testing[12]=1
            }
          }
          else{
            testing[13]=1
          }
        }
        else{
          testing[14]=1
        }
      }
      
      # find out which interaction hypothesis of 2 factors is tested
      else if(l==2){
        if(grepl(group,x[i])){
          if(grepl(subgroup,x[i])){
            testing[5]=1
          }
          if(grepl(factor1,x[i])){
            testing[6]=1
          }
          if(grepl(factor2,x[i])){
            testing[7]=1
          }
        }
        else if(grepl(subgroup,x[i])){
          if(grepl(factor1,x[i])){
            testing[8]=1
          }
          if(grepl(factor2,x[i])){
            testing[9]=1
          }
        }
        else if(grepl(factor1,x[i])){
          if(grepl(factor2,x[i])){
            testing[10]=1
          }
        }
      }
      # l = 1
      else{
        if(group == x[i]){testing[1]=1}
        else if(subgroup == x[i]){testing[2]=1}
        else if(factor1 == x[i]){testing[3]=1}
        else {testing[4]=1}
      }
    }
    X=data
    data = colnames(dat)[1]
    return(hrm.test.2.between.within(X, alpha, group , subgroup, factor1, factor2, subject, data, testing ))
  }
}