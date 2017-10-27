####################################################################################################################################
### Filename:    Repeated.R
### Description: Main file, where end-user functions are defined;
###              EEG data example
###
###
####################################################################################################################################

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
hrm.test <- function(n, a, d, data, alpha=0.05){
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
hrm.test.matrix <- function(data, alpha=0.05){
  
  if(!is.list(data)){
    stop("data needs to be a list containing the data matrices of all groups")
  }
  
  a <- length(data)
  if(a < 2){
    stop("At least two groups are needed.")
  }
  
  n <- rep(0,a)
  for(i in 1:a){
    tmp <- data[[i]]
    if(!is.matrix(tmp)){
      stop("The elements of data need to be matrices.")
    }
    n[i] <- dim(tmp)[1]
  }
  d <- dim(data[[1]])[2]
  d2 <- dim(data[[2]])[2]
  if(d != d2){
    stop("The number of measurements for each group need to be the same.")
  }
  
  stopifnot(is.list(data), all(is.finite(n)), alpha<=1, alpha>=0, d>=1, a>=1)
  
  temp0 <- hrm.A.weighted(n,a,d,data,alpha)
  temp1 <- hrm.A.unweighted(n,a,d,data,alpha)
  temp2 <- hrm.B(n,a,d,data,alpha)
  temp3 <- hrm.AB(n,a,d,data,alpha)
  temp4 <- hrm.A_B(n,a,d,data,alpha)

  output <- list()
  output$result <- rbind(temp0,temp1,temp2,temp3,temp4)
  output$formula <- NULL
  output$alpha <- alpha
  output$subject <- NULL
  output$factors <- list(NULL, NULL)
  class(output) <- "HRM"
  
  return (output)
}

#' Test for no main effects and interactino effects of one between-subject factor and one crossed within-subject factors
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.one <- function(X, alpha, group , factor1, subject, data, testing = rep(1,4), formula ){

  temp0 <- if(testing[1]) {hrm.1w.1f(X, alpha, group , factor1,  subject, data, "A", paste(as.character(group), " (weighted)"))}
  temp1 <- if(testing[2]) {hrm.1w.1f(X, alpha, group , factor1,  subject, data, "Au", paste(as.character(group)))}
  temp2 <- if(testing[3]) {hrm.1w.1f(X, alpha, group , factor1,  subject, data, "B", paste(as.character(factor1)))}
  temp3 <- if(testing[4]) {hrm.1w.1f(X, alpha, group , factor1, subject, data, "AB", paste(as.character(group), ":",as.character(factor1)))}
  
  output <- list()
  output$result <- rbind(temp0, temp1, temp2, temp3)
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group), c(factor1))
  class(output) <- "HRM"
  
  return (output)
}

#' Test for no main effects and interactino effects of one between-subject factor and two crossed within-subject factors
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param factor2 column name of the data frame X of the second factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.within <- function(X, alpha, group , factor1, factor2, subject, data, testing = rep(1,7), formula ){
  
  # create list for storing results; NULL used, because it is ignored by rbind
  temp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  for(i in 1:7){
    if(testing[i]) {
      temp[[i]] <- hrm.1w.2f(X, alpha, group , factor1, factor2, subject, data, H = i )
    }
  }
  
  output <- list()
  output$result <- rbind(temp[[1]], temp[[2]], temp[[3]], temp[[4]], temp[[5]], temp[[6]], temp[[7]])
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group), c(factor1, factor2))
  class(output) <- "HRM"
  
  return (output)
}



#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor column name of the data frame X of within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.between <- function(X, alpha, group , subgroup, factor, subject, data, testing = rep(1,7), formula ){

  # create list for storing results; NULL used, because it is ignored by rbind
  temp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  for(i in 1:7){
    if(testing[i]) {
      temp[[i]] <- hrm.2w.1f(X, alpha, group , subgroup, factor, subject, data, H = i )
    }
  }
  
  output <- list()
  output$result <- rbind(temp[[1]], temp[[2]], temp[[3]], temp[[4]], temp[[5]], temp[[6]], temp[[7]])
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group, subgroup), c(factor))
  class(output) <- "HRM"
 
  return (output)
}


#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param subgroup column name of the subgroups (crossed with groups)
#' @param factor1 column name of the data frame X of the first within-subject factor
#' @param factor2 column name of the data frame X of the second within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.2.between.within <- function(X, alpha, group , subgroup, factor1, factor2, subject, data, testing = rep(1,15), formula ){
  
  # create list for storing results; NULL used, because it is ignored by rbind
  temp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  for(i in 1:15){
    if(testing[i]) {
      temp[[i]] <- hrm.2w.2f(X, alpha, group , subgroup, factor1, factor2, subject, data, H = i )
    }
  }
  output <- list()
  output$result <- rbind(temp[[1]], temp[[2]], temp[[3]], temp[[4]], temp[[5]], temp[[6]], temp[[7]], temp[[8]], temp[[9]], temp[[10]], temp[[11]], temp[[12]], temp[[13]], temp[[14]], temp[[15]])
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group, subgroup), c(factor1, factor2))
  class(output) <- "HRM"
  return (output)
}


#' Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first within-subject factor
#' @param factor2 column name of the data frame X of the second within-subject factor
#' @param factor3 column name of the data frame X of the third within-subject factor
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the data frame X containing the measurement data
#' @param testing vector specifying which hypotheses should be tested
#' @param formula formula object from the user input
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.test.3.between <- function(X, alpha, group , factor1, factor2, factor3, subject, data, testing = rep(1,15), formula ){
  temp0 <- if(testing[1]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3, subject, data, P, J, J, J,  paste(as.character(group) ) )}
  temp1 <- if(testing[2]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, P, J, J, paste(as.character(factor1)) )}
  temp2 <- if(testing[3]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, J, P, J, paste(as.character(factor2)) )}
  temp3 <- if(testing[4]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, J, J, P, paste(as.character(factor3)) )}
  temp4 <- if(testing[5]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, P, J, J, paste(as.character(group),":",as.character(factor1)) )} 
  temp5 <- if(testing[6]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, J, P, J, paste(as.character(group),":",as.character(factor2)) )}
  temp6 <- if(testing[7]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, J, J, P, paste(as.character(group),":",as.character(factor3)) )}
  temp7 <- if(testing[8]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, P, P, J, paste(as.character(factor1),":",as.character(factor2)) )}
  temp8 <- if(testing[9]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, P, J, P, paste(as.character(factor1),":",as.character(factor3)) )}
  temp9 <- if(testing[10]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, J, P, P,paste(as.character(factor2),":",as.character(factor3)) )}
  temp10 <- if(testing[11]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, P, P, J, paste(as.character(group),":",as.character(factor1), ":", as.character(factor2)) )}
  temp11 <- if(testing[12]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, P, J, P, paste(as.character(group),":",as.character(factor1), ":", as.character(factor3)) )}
  temp12 <- if(testing[13]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, J, P, P, paste(as.character(group),":",as.character(factor2), ":", as.character(factor3)) )}
  temp13 <- if(testing[14]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, J, P, P, P, paste(as.character(factor1),":",as.character(factor2), ":", as.character(factor3)) )}
  temp14 <- if(testing[15]) {hrm.1w.3f(X, alpha, group, factor1, factor2, factor3,  subject, data, P, P, P, P, paste(as.character(group),":",as.character(factor1), ":", as.character(factor2), ":", as.character(factor3)) )}
  
  output <- list()
  output$result <- rbind(temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14)
  output$formula <- formula
  output$alpha <- alpha
  output$subject <- subject
  output$factors <- list(c(group), c(factor1, factor2, factor3))
  class(output) <- "HRM"
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
hrm.test.2 <- function(X, alpha = 0.05, group , subgroup, factor1, factor2, subject, data ){
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
#' @param factor2 column name within the data frame data specifying the the second subplot-factor
#' @param factor3 column name within the data frame data specifying the the third subplot-factor
#' @param subject column name within the data frame X identifying the subjects
#' @param response column name within the data frame X containing the response variable
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @example R/example_2.txt
#' @keywords export
hrm.test.dataframe <- function(data, alpha = 0.05, group , subgroup, factor1, factor2, factor3, subject, response ){
  
  temp <- data
  data <- response
  X <- temp
  
  # if(missing(factor1) || !is.character(factor1)){
  #   print("At least one within-subject factor is needed!")
  #   stop("factor1 column name not specified ")
  # }
  
  # if(missing(group) || !is.character(group)){
  #   print("At least one between-subject factor is needed!")
  #   stop("group column name not specified ")
  # }
  
  if(!missing(factor1) & !missing(factor2) & !missing(factor3) & !missing(group)  & !missing(subgroup)  ){
    stop("The maximum number of factors that can be used is four.")
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
  dat <- X
  dat <- data.frame(dat, subj = X[,subject])
  s1<-subset(dat, dat$subj==dat$subj[1])
  measurements <- dim(s1)[1]
  countSubplotFactor <- nlevels(s1[,factor1])
  if(!missing(factor2)){
    countSubplotFactor <- countSubplotFactor*nlevels(s1[,factor2])
  }
  if(!missing(factor3)){
    countSubplotFactor <- countSubplotFactor*nlevels(s1[,factor3])
  }
  if(countSubplotFactor==0 & measurements == 1){
    countSubplotFactor <- 1
  }
  if(!(measurements == countSubplotFactor)){
    stop(paste("The number of repeated measurements per subject (", measurements, ") is uneqal to the number of levels of the subplot factors (", countSubplotFactor, ")."))
  }
  
  # case: no subplot, one wholeplot
  if(missing(factor2) & missing(subgroup) & missing(factor3) & missing(factor1) & !missing(group)){
    
    if(!is.factor(X[,group])){
      stop(paste("The column ", group, " is not a factor." ))
    }
    
    return(hrm.test.1.none(X, alpha , group, subject, data, formula = NULL ))
  }

  # case: one subplot, no wholeplot
  if(missing(factor2) & missing(subgroup) & missing(factor3) & missing(group) & !missing(factor1)){
    
    if(!is.factor(X[,factor1])){
      stop(paste("The column ", factor1, " is not a factor." ))
    }
    
    return(hrm.test.1.one(X, alpha , factor1, subject, data, formula = NULL ))
  }

  
  if(missing(factor2) & missing(subgroup) & missing(factor3)){
    if(!is.factor(X[,group])){
      stop(paste("The column ", group, " is not a factor." ))
    }
    if(!is.factor(X[,factor1])){
      stop(paste("The column ", factor1, " is not a factor." ))
    }
    
    return(hrm.test.2.one(X, alpha, group , factor1, subject, data, formula = NULL ))
  }
  
  
  if(missing(factor2) & !missing(subgroup) & missing(factor3)){
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
    return(hrm.test.2.between(X, alpha, group , subgroup, factor1, subject, data, formula = NULL))
  }
  
  if(!missing(factor2) & missing(subgroup) & missing(factor3)){
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
    return(hrm.test.2.within(X, alpha, group , factor1, factor2, subject, data, formula = NULL ))
  }
  
  if(!missing(factor2) & !missing(subgroup) & missing(factor3)){
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
    return(hrm.test.2.between.within(X, alpha, group , subgroup, factor1, factor2, subject, data, formula = NULL))
  }
  if(!missing(factor1) & !missing(factor2) & !missing(factor3) & !missing(group)  & missing(subgroup)  ){
    if(!is.character(factor1)){
      stop("factor1 column name not specified")
    }
    if(!is.character(factor2)){
      stop("factor2 column name not specified")
    }
    if(!is.character(factor3)){
      stop("factor3 column name not specified")
    }
    if(!is.character(group)){
      stop("group column name not specified")
    }
    if(!is.factor(X[,group])){
      stop(paste("The column ", group, " is not a factor." ))
    }
    if(!is.factor(X[,factor1])){
      stop(paste("The column ", factor1, " is not a factor." ))
    }
    if(!is.factor(X[,factor2])){
      stop(paste("The column ", factor2, " is not a factor." ))
    }
    if(!is.factor(X[,factor3])){
      stop(paste("The column ", factor3, " is not a factor." ))
    }
    formula <- as.formula(paste(data, "~", group, "*", factor1, "*", factor2, "*", factor3))
    
    return(hrm_test(formula=formula,alpha=alpha,subject=subject, data=X ))
    
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
hrm_test <- function(formula, data, alpha = 0.05,  subject ){
  
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
  dat <- model.frame(formula, data)
  dat2 <- data.frame(dat,subj=data[,subject])
  m <- ncol(dat)
  
  if(!is.numeric(dat[,1])){
    stop("Response variable needs to be numeric!")
  }
  
  
  
  # find out, in which columns are the wholeplot or subplot factors
  s1<-subset(dat2, dat2$subj==dat2$subj[1])
  measurements <- dim(s1)[1]
  countSubplotFactor <- 1
  wholeplot<-rep(-1, m)
  subplot<-rep(-1, m)
  for(i in 2:m){
    if(!is.factor(dat2[,i])){
      stop(paste("The column ", colnames(dat2)[i], " is not a factor." ))
    }
    if(length(unique(s1[,i]))==nlevels(dat2[,i])){
      subplot[i]<-1
      countSubplotFactor <- countSubplotFactor*nlevels(s1[,i])
    }
    else{
      wholeplot[i]<-1
    }
  }
  wholeplot <- which(wholeplot==1)
  subplot <- which( subplot==1)

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
  if(length(subplot)<1 & length(wholeplot)>1){
    stop("The model needs at least one within-subject factor.")
  }
  
  # Case: no wholeplot, one subpot factor
  if(length(wholeplot) < 1 & length(subplot) == 1){
    factor1 <- colnames(dat2)[subplot[1]]
    x<-attributes(terms.formula(formula))$term.labels
    
    subplot<-colnames(dat2[,subplot])
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.1.one(X, alpha , factor1, subject, data, formula ))   
    
  }
  
  # Case: one wholeplot, no subpot factor
  if(length(wholeplot) == 1 & length(subplot) < 1){
    group <- colnames(dat2)[wholeplot[1]]
    x<-attributes(terms.formula(formula))$term.labels
    
    wholeplot<-colnames(dat2[,group])
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.1.none(X, alpha , group, subject, data, formula ))  
  
  }
  
  # Case: 1 whole and 1 subplot factor
  if(length(wholeplot)==1 & length(subplot)==1){
    group <- colnames(dat2)[wholeplot[1]]
    factor1 <- colnames(dat2)[subplot[1]]
    x<-attributes(terms.formula(formula))$term.labels
    
    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])
    
    testing <- rep(0,4)
    for(i in 1:length(x)){
      
      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])
      
      # interaction hypothesis of 2 factors
      if(l == 2){testing[4]<-1}
      
      else{
        if(group == x[i]){
          testing[1]<-1
          testing[2]<-1
        }
        else if(factor1 == x[i]){testing[3]<-1}
      }
    }
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.2.one(X, alpha, group , factor1, subject, data, testing, formula ))   
  }
  
  # Case: 2 wholeplot, 1 subplot factor
  if(length(wholeplot)==2 & length(subplot)==1){
    group <- colnames(dat2)[wholeplot[1]]
    subgroup <- colnames(dat2)[wholeplot[2]]
    factor1 <- colnames(dat2)[subplot[1]]
    x<-attributes(terms.formula(formula))$term.labels
    
    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])
    
    testing <- rep(0,7)
    for(i in 1:length(x)){
      
      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])
      
      # interaction hypothesis of 4 factors
      if(l == 3){testing[7]<-1}
      
      # find out which interaction hypothesis of 2 factors is tested
      else if(l==2){
        if(grepl(group,x[i])){
          if(grepl(subgroup,x[i])){
            testing[4]<-1
          }
          if(grepl(factor1,x[i])){
            testing[5]<-1
          }
        }
        else if(grepl(subgroup,x[i])){
          if(grepl(factor1,x[i])){
            testing[6]<-1
          }
        }
      }
      # l = 1
      else{
        if(group == x[i]){testing[1]<-1}
        else if(subgroup == x[i]){testing[2]<-1}
        else if(factor1 == x[i]){testing[3]<-1}
      }
    }
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.2.between(X, alpha, group , subgroup, factor1, subject, data, testing, formula )) 
  }
  
  # Case: 1 wholeplot, 2 subplot factors
  if(length(wholeplot)==1 & length(subplot)==2){
    group <- colnames(dat2)[wholeplot[1]]
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    x<-attributes(terms.formula(formula))$term.labels
    
    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])
    
    testing <- rep(0,7)
    for(i in 1:length(x)){
      
      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])
      
      # interaction hypothesis of 3 factors
      if(l == 3){testing[7]<-1}
      
      # find out which interaction hypothesis of 2 factors is tested
      else if(l==2){
        if(grepl(group,x[i])){
          if(grepl(factor1,x[i])){
            testing[4]<-1
          }
          if(grepl(factor2,x[i])){
            testing[5]<-1
          }
        }
        else if(grepl(factor1,x[i])){
          testing[6]<-1
        }
      }
      # l = 1
      else{
        if(group == x[i]){testing[1]<-1}
        else if(factor1 == x[i]){testing[2]<-1}
        else {testing[3]<-1}
      }
    }
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.2.within(X, alpha, group , factor1, factor2, subject, data, testing, formula ))
  }
  
  # Case: 1 wholeplot, 3 subplot factors
  if(length(wholeplot)==1 & length(subplot)==3){
    group <- colnames(dat2)[wholeplot[1]]
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    factor3 <- colnames(dat2)[subplot[3]]
    
    x<-attributes(terms.formula(formula))$term.labels
    
    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])
    
    testing <- rep(0,15)
    for(i in 1:length(x)){
      
      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])
      
      # interaction hypothesis of 3 factors
      if(l == 4){testing[15]<-1}
      
      # find out which interaction hypothesis of 3 factors is tested
      else if(l==3){
        if(grepl(group,x[i])){
          if(grepl(factor1,x[i])){
            if(grepl(factor2,x[i])){
              testing[11]<-1
            }
            else{
              testing[12]<-1
            }
          }
          if(grepl(factor2,x[i])){
            if(grepl(factor3,x[i])){
              testing[13]<-1
            }
          }
        }
        else if(grepl(factor1,x[i])){
          testing[14]<-1
        }
      }
      # l = 2
      else if (l==2){
        if(grepl(group,x[i])){
          if(grepl(factor1,x[i])){
            testing[5]<-1
          }
          if(grepl(factor2,x[i])){
            testing[6]<-1
          }
          if(grepl(factor3,x[i])){
            testing[7]<-1
          }
        }
        else if(grepl(factor1,x[i])){
          if(grepl(factor2,x[i])){
            testing[8]<-1
          }
          if(grepl(factor3,x[i])){
            testing[9]<-1
          }
        }
        else if(grepl(factor2,x[i])){
          if(grepl(factor3,x[i])){
            testing[10]<-1
          }
        }
      }
      else if(l==1){
        if(grepl(group, x[[i]])){
          testing[1]<-1
        }
        if(grepl(factor1, x[[i]])){
          testing[2]<-1
        }
        if(grepl(factor2, x[[i]])){
          testing[3]<-1
        }
        if(grepl(factor3, x[[i]])){
          testing[4]<-1
        }
      }
    }
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.3.between(X, alpha, group , factor1, factor2, factor3, subject, data, testing, formula ))
  }
  
  
  
  if(length(wholeplot)==2 & length(subplot)==2){
    group <- colnames(dat2)[wholeplot[1]]
    subgroup <- colnames(dat2)[wholeplot[2]]
    factor1 <- colnames(dat2)[subplot[1]]
    factor2 <- colnames(dat2)[subplot[2]]
    x<-attributes(terms.formula(formula))$term.labels
    
    wholeplot<-colnames(dat2[,wholeplot])
    subplot<-colnames(dat2[,subplot])
    
    testing <- rep(0,15)
    for(i in 1:length(x)){
      
      tmp <- strsplit(x[i],":")
      l <- length(tmp[[1]])
      
      # interaction hypothesis of 4 factors
      if(l == 4){testing[15]<-1}
      
      # find out which interaction hypothesis of 3 factors is tested
      else if(l==3){
        if(grepl(group,x[i])){
          if(grepl(subgroup,x[i])){
            if(grepl(factor1,x[i])){
              testing[11]<-1
            }
            else{
              testing[12]<-1
            }
          }
          else{
            testing[13]<-1
          }
        }
        else{
          testing[14]<-1
        }
      }
      
      # find out which interaction hypothesis of 2 factors is tested
      else if(l==2){
        if(grepl(group,x[i])){
          if(grepl(subgroup,x[i])){
            testing[5]<-1
          }
          if(grepl(factor1,x[i])){
            testing[6]<-1
          }
          if(grepl(factor2,x[i])){
            testing[7]<-1
          }
        }
        else if(grepl(subgroup,x[i])){
          if(grepl(factor1,x[i])){
            testing[8]<-1
          }
          if(grepl(factor2,x[i])){
            testing[9]<-1
          }
        }
        else if(grepl(factor1,x[i])){
          if(grepl(factor2,x[i])){
            testing[10]<-1
          }
        }
      }
      # l = 1
      else{
        if(group == x[i]){testing[1]<-1}
        else if(subgroup == x[i]){testing[2]<-1}
        else if(factor1 == x[i]){testing[3]<-1}
        else {testing[4]<-1}
      }
    }
    X<-data
    data <- colnames(dat)[1]
    return(hrm.test.2.between.within(X, alpha, group , subgroup, factor1, factor2, subject, data, testing, formula ))
  }
}