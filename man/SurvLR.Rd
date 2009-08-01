\name{SurvLR}
\alias{SurvLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transform Surv object to data frame with L and R values}
\description{
Takes a \code{\link[survival]{Surv}} object and transforms it into a data frame with
two variables, L and R, representing the left and right interval of interval censored 
data. The failure time is known to be in the interval (L,R]. Right censored data 
are handled by setting L=R for observed and R=Inf for right censored.  These are interpreted correctly 
by \code{\link{icfit}} and \code{\link{ictest}}.
}
\usage{
SurvLR(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{\link[survival]{Surv}} object}
}
\details{Currently type='counting' not supported.}
\value{
  A data frame with two variables:
  \item{L}{left end of interval}
  \item{R}{right end of interval}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
\seealso{ Called from \code{\link{icfit}} and \code{\link{ictest}} }
\examples{
time<-c(1,5,3,7)
status<-c(1,1,0,1)
y<-Surv(time,status)
SurvLR(y)
}
\keyword{misc}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
