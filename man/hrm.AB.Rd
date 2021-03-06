% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f2_old.R
\name{hrm.AB}
\alias{hrm.AB}
\title{Test for no interaction between treatment and time}
\usage{
hrm.AB(n, a, d, X, alpha, nonparametric = FALSE)
}
\arguments{
\item{n}{an vector containing the sample sizes of all groups}

\item{a}{number of groups}

\item{d}{number of dimensions (time points)}

\item{X}{list containing the data matrices of all groups}

\item{alpha}{alpha level used for the test}
}
\value{
Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
}
\description{
Test for no interaction between treatment and time
}
\keyword{internal}
