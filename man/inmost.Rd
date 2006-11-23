\name{inmost}
\alias{inmost}
\title{Innermost Intervals for Interval-Censored Data}
\description{
This function calculates the innermost intervals (Turnbull's
algorithm) for interval-censored data.  Right-censored data
is not allowed at this point.
}
\usage{
inmost(data, eps)
}
\arguments{
  \item{data}{A matrix with two columns, consisting of left and right
endpoints of the interval data}
  \item{eps}{A scalar additive value which jitters the data slightly to
break ties.}
}
\value{A 2 column matrix of the left and right endpoints of the
innermost intervals.}
\references{Turnbull (1976) The empirical distribution function 
with arbitrarily grouped, censored and truncated data. Journal
of the Royal Statistical Society Series B, 38: 290-295.}
\author{W.J. Braun}
\keyword{models}
