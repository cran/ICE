\name{likelihoodcv}
\alias{likelihoodcv}
\title{Log Cross-Validated Likelihood}
\description{
Log Cross-Validated Likelihood for interval-censored data.  The
likelihood is the product of integrals over the innermost intervals.
Leave-one-out cross-validation here is accomplished by leaving
out each innermost interval and re-computing the integral using the
remaining data.
}
\usage{
likelihoodcv(p, data, m=101)
}
\arguments{
  \item{p}{A scalar specifying the bandwidth used in the kernel density
estimate.}
  \item{data}{A matrix with two columns, consisting of left and right
endpoints of the interval data}
  \item{m}{numeric quantity which controls the number of grid points
used in the density estimate}
}
\value{A scalar value.}
\references{Braun, J., Duchesne, T. and Stafford, J.E. (2005)
Local likelihood density estimation for interval censored
data.  Canadian Journal of Statistics 33: 39-60.}
\author{W.J. Braun}
\keyword{models}
