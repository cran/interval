\name{bcos}
\alias{bcos}
\docType{data}
\title{Breast Cosmesis Data}
\description{
The often used data set for interval censored data, described and given in full in Finkelstein and Wolfe (1985).
}
\usage{data(bcos)}
\format{
  A data frame with 94 observations on the following 3 variables.
  \describe{
    \item{\code{left}}{a numeric vector}
    \item{\code{right}}{a numeric vector}
    \item{\code{treatment}}{a factor with levels \code{Rad} and \code{RadChem}}
  }
}
%\details{}
\source{
Finkelstein, D.M., and Wolfe, R.A. (1985). A semiparametric model for regression analysis of interval-censored 
failure time data. Biometrics 41: 731-740.
}
%\references{ ~~ possibly secondary sources and usages ~~}
\examples{
data(bcos)
}
\keyword{datasets}
