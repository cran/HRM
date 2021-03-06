% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Repeated.R
\name{hrm.test.2.between.within}
\alias{hrm.test.2.between.within}
\title{Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor}
\usage{
hrm.test.2.between.within(
  X,
  alpha,
  group,
  subgroup,
  factor1,
  factor2,
  subject,
  data,
  testing = rep(1, 15),
  formula,
  nonparametric,
  np.correction
)
}
\arguments{
\item{X}{dataframe containing the data in the long table format}

\item{alpha}{alpha level used for the test}

\item{group}{column name of the data frame X specifying the groups}

\item{subgroup}{column name of the subgroups (crossed with groups)}

\item{factor1}{column name of the data frame X of the first within-subject factor}

\item{factor2}{column name of the data frame X of the second within-subject factor}

\item{subject}{column name of the data frame X identifying the subjects}

\item{data}{column name of the data frame X containing the measurement data}

\item{testing}{vector specifying which hypotheses should be tested}

\item{formula}{formula object from the user input}
}
\value{
Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
}
\description{
Test for no main effects and interaction effects of two crossed between-subject factors and one within-subject factor
}
\keyword{internal}
