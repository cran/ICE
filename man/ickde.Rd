\name{ickde}
\alias{ickde}
\title{Interval-Censored Kernel Density Estimation}
\description{
Iterated conditional expectation kernel density estimation 
using a local constant.  The bandwidth is assumed fixed. (See
the example for a way to get a quick ballpark estimate of the bandwidth.)
The gaussian, epanechnikov and biweight kernels can be used. Note
that the bandwidth estimate would have to be adjusted before using
with epanechnikov or biweight. 
}
\usage{
ickde(I, h, f, m, n.iterations = 10, x1, xm, right.limit = 10000,kernel="gaussian")
}
\arguments{
  \item{I}{A matrix with two columns, consisting of left and right
endpoints of the interval data}
  \item{h}{A scalar bandwidth}
  \item{f}{An initial estimate of the density at a sequence of grid 
points (optional; if this is used, do not specify m)}
  \item{m}{The number of (equally-spaced) grid points at which
the density is to be estimated}
  \item{n.iterations}{The maximum number of iterations allowed}
  \item{x1}{The left-most grid point (optional)}
  \item{xm}{The right-most grid point (optional)}
  \item{right.limit}{For right-censored data, the value to be
used as an artificial right endpoint for the intervals}
  \item{kernel}{character argument indicated choice of kernel; current
choices are "gaussian", "epanechnikov", "biweight"}
}
\value{An object of class \code{IC}
}
\references{Braun, J., Duchesne, T. and Stafford, J.E. (2005)
Local likelihood density estimation for interval censored
data.  Canadian Journal of Statistics 33: 39-60.}
\author{W.J. Braun}
\seealso{\code{\link{dpik}}}
\examples{
 tmp <- apply(ICHemophiliac, 1, mean)
 h <- try(dpik(tmp), silent=T) # dpik() will work if KernSmooth is loaded
 if (class(h) !="numeric" ) h <- .9  # this makes the example work 
                       # if KernSmooth is not loaded
 estimate <- ickde(ICHemophiliac, m=200, h=h)
 plot(estimate, type="l")
}
\keyword{models}
