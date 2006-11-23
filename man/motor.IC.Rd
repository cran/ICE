\name{motor.IC}
\alias{motor.IC}
\docType{data}
\title{Motorcycle Data with artificial interval-censored responses}
\description{
The original motorcycle data with exponential random variables
added and subtracted to the responses to create interval-censored
responses.
}
\usage{data(motor.IC)}
\format{
  A data frame with 4 variables
  \describe{
    \item{V1}{time, the explanatory vector}
    \item{V2}{acceleration, the response}
    \item{y.L}{y - an exponential random variable}
    \item{y.R}{y + an exponential random variable}
  }
}
\keyword{datasets}
