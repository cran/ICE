\name{iclocpoly}
\alias{iclocpoly}
\title{Interval-Censored Local Polynomial Regression Estimation}
\description{Local polynomial regression estimation for
interval-censored data.
}
\usage{
iclocpoly(x, y=NULL, y.IC, degree=0, h, niter=10, kernel="normal", gridsize=401)
}
\arguments{
\item{x}{uncensored explanatory variable vector}
\item{y}{uncensored portion of response vector (optional)}
\item{y.IC}{two-column matrix of left and right interval
endpoints for censored responses}
\item{degree}{degree of local polynomial}
\item{h}{bandwidth}
\item{niter}{number of iterations}
\item{kernel}{smoothing kernel to be used; default is "normal"; other
choices as in the function locpoly}
\item{gridsize}{number of gridpoints; again as in locpoly}
}
\value{A list consisting
of the explanatory variable x and the imputed responses y as well
as the estimate of sigma.  Usually, one would apply locpoly with
a similar bandwidth to obtain the final fitted model.
}
\examples{
library(KernSmooth)
data(motor.IC)
estimate <- iclocpoly(x=motor.IC$V1, y.IC = cbind(motor.IC$y.L, 
motor.IC$y.R), degree=0, h=1)
plot(motor.IC$V1, motor.IC$V2)
lines(locpoly(estimate$x, estimate$y, bandwidth=1, degree=0),col=4)
}
\keyword{models}
