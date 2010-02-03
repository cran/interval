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
right: take the right side  of the non-unique MLE interval. The LRin attribute is ignored (see warning).

If \code{icfit} has more than one strata, then performs the operations on each 
stratum.
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
if there is only one stratum, then creates a LIST,with elements
  \item{S}{vector of survival function values at each element of times}
  \item{times}{vector of times for which need survival function}
  \item{unique.mle}{logical denoting whether associated survival value is a unique MLE}
  \item{nonUMLE.method}{character vector describing non-unique MLE method}
if there are k>1 strata, then creates a list with k+1 elements, the elements 1:k are 
lists of results for each strata, and element k+1 is called strataNames and is a character 
vector of strata names.  
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\section{Warning}{The \code{getsurv} function does not use LRin attributes, so values exactly on the 
intmap values may only represent the limit approaching that value, not the survival at that value.}
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
