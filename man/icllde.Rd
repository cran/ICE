\name{icllde}
\alias{icllde}
\title{Interval-Censored Local Linear Density Estimation}
\description{
This is the local linear version of \code{ickde}.
}
\usage{
icllde(I, h, f, m, n.iterations = 10, x1, xm, right.limit = 10000, kernel="gaussian")
}
\arguments{
  \item{I}{A matrix with two columns, consisting of left and right
          endpoints of the interval data}
  \item{h}{A scalar bandwidth}
  \item{f}{An initial estimate of the density at a sequence of grid 
          points (optional)}
  \item{m}{The number of (equally-spaced) grid points at which the
          density is to be estimated}
  \item{n.iterations}{The maximum number of iterations allowed}
  \item{x1}{The minimum grid point (optional)}
  \item{xm}{The maximum grid point (optional)}
  \item{right.limit}{For right-censored data, the value to be used as an
          artificial right endpoint for the intervals}
  \item{kernel}{character argument indicated choice of kernel; current
choices are "gaussian", "epanechnikov", "biweight"}
}
\value{An object of class \code{IC}
 }
\examples{
     data(ICHemophiliac)
     estimate <- icllde(ICHemophiliac, m=200, h=.9)
     plot(estimate, type="l")
}
\keyword{models}
