\name{Aintmap}
\alias{Aintmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create A matrix and intmap}
\description{
The A matrix is an n by k matrix of zeros and ones, where each row represents one of n failure times, 
and each column represents 
a possible interval for the nonparametric maximum likelihood estimate (NPMLE). 
The function \code{Aintmap} creates an A matrix and associated intmap 
from left and right intervals (L and R) which may not may not include the boundary 
of the interval (using Lin or Rin).
The matrix intmap  denotes the intervals of the potential jumps in the distribution of the 
NPMLE, 
and its attribute LRin denotes whether to include each of the intervals 
or not. Called  by \code{\link{icfit}}.
}
\usage{
Aintmap(L,R,Lin=NULL,Rin=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{numeric vector of left endpoints of censoring interval}
  \item{R}{numeric vector of right endpoints of censoring interval}
  \item{Lin}{logical vector, should L be included in the interval? (see details)}
  \item{Rin}{logical vector, should R be included in the interval? (see details)} 
}
%\details{}
\value{
  A list with two objects:
  \item{A}{an n by k matrix of 0 and 1s}
  \item{intmap}{the associated intmap}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
\seealso{ Called from \code{\link{icfit}} and \code{\link{ictest}} }
\examples{
Aintmap(c(2,3,3,7),c(3,5,5,8),Lin=c(FALSE,TRUE,FALSE,FALSE),Rin=c(TRUE,FALSE,TRUE,FALSE))
}
\keyword{misc}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
