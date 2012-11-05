\name{bickde}
\alias{bickde}
\title{Bandwidth choice for Interval-Censored Kernel Density Estimation}
\description{
Likelihood Cross-Validation bandwidth choice for interval-censored
kernel density estimates.  Also computed is the direct-plug-in
estimate (using the KernSmooth function dpik based
on the interval midpoints.
}
\usage{
bickde(data, factor)
}
\arguments{
  \item{data}{A matrix with two columns, consisting of left and right
endpoints of the interval data}
  \item{factor}{A scalar factor which gives upper and lower
bounds for the initial interval for the golden section search relative
to a preliminary value computed by dpik.}
Default value is 10 which specifies an initial interval of (dpik/10, dpik*10).}
\value{A scalar value.}
\details{Maximization of the likelihood is accomplished by
golden section search using the optimize() function.  This
can be very SLOW.}
\references{Braun, J., Duchesne, T. and Stafford, J.E. (2005)
Local likelihood density estimation for interval censored
data.  Canadian Journal of Statistics 33: 39-60.}
\author{W.J. Braun}
\examples{
# Not run
# bickde(ICHemophiliac)
}
\keyword{models}
