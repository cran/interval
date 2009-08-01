\name{getsurv}
\alias{getsurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get survival values from  icfit object}
\description{
For a vector of times, \code{getsurv} gets the associated survival values. 
 The MLE is not uniquely defined for times inbetween the first and second row on the 
 same column of the intmap. 
If there is not 
a unique MLE for a specific time, then either use, 
interpolation: (default), which basically finds the point on the line connecting 
the two points bounding the non-unique MLE interval, or,
left: take the left side of the non-unique MLE interval (smallest value) or,
right: take the right side  of the non-unique MLE interval. 
}
\usage{
getsurv(times, icfit, nonUMLE.method = "interpolation")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{numeric vector of times }
  \item{icfit}{ icfit object used to define the survival function}
  \item{nonUMLE.method}{character vector, either "interpolation", "left" or "right".
Method for finding survival when times element is not 
at a unique MLE time.}
}
%\details{}
\value{
a LIST,with elements
  \item{S}{vector of survival function values at each element of times}
  \item{times}{vector of times for which need survival function}
  \item{unique.mle}{logical denoting whether associated survival value is a unique MLE}
  \item{nonUMLE.method}{character vector describing non-unique MLE method}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
