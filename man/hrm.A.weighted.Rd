% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f2_old.R
\name{hrm.A.weighted}
\alias{hrm.A.weighted}
\title{Test for no main treatment effect (weighted version)}
\usage{
hrm.A.weighted(n, a, d, X, alpha, nonparametric = FALSE)
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
Test for no main treatment effect (weighted version)
}
\keyword{internal}
