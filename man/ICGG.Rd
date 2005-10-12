\name{ICGG}
\alias{ICGG}
\docType{data}
\title{Gentleman and Geyer's Data}
\description{
Gentleman and Geyer's Data.
}
\usage{data(ICGG)}
\format{
  A data frame with 6 observations on the following 2 variables.
  \describe{
    \item{left}{the left interval endpoint}
    \item{right}{the right interval endpoint}
  }
}
\source{
Gentleman, R. and Geyer, C. (1994) Maximum Likelihood for Interval
Censored Data: Consistency and Computation.  Biometrika. 
}
\examples{
 library(KernSmooth)
 tmp <- apply(ICGG, 1, mean)  # tmp now contains the interval midpoints
 h <- dpik(tmp)               # direct-plug-in bandwidth selected for
                              # interval midpoints
 par(mfrow=c(1,2))
 estimate <- ickde(ICGG, h=h, m=200)
 plot(estimate, type="l", main="One fixed point")
 estimate <- ickde(ICGG, f=c(rep(1,60),rep(0,90),rep(1,60)), h=.1)
 plot(estimate, type="l", main="Another fixed point")
}
\keyword{datasets}
