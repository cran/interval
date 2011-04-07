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
\details{
The \code{Lin} and \code{Rin} specify whether or not to include the ends of the intervals.
They may be length 1 (and apply to all n values) or length n. The function automatically 
only returns the innermost intervals (also called the Turnbull intervals [see Turnbull, 1976], or the regions of the maximal cliques
[see Gentleman and Vandal, 2002]).
The innermost intervals give the "primary reduction" of Aragon and Eberly (1992). 
}
\value{
  A list with two objects:
  \item{A}{an n by k matrix of 0 and 1s}
  \item{intmap}{the associated intmap}
}
\references{ 
   Aragon, J and Eberly, D (1992). On convergence of convex 
    minorant algorithms for distribution estimation with 
    interval-censored data.  J. of Computational and Graphical 
    Statistics. 1: 129-140.

Gentleman R, and Vandal, A (2002). Nonparametric estimation of the bivariate CDF for arbitrarily censored dtaa. Canadian J of Stat 30: 557-571.

 Turnbull, B (1976). The empirical distribution function with arbitrarily grouped, censored and truncated data.
JRSS-B, 38: 290-295.


}
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
\seealso{ Called from \code{\link{icfit}} and \code{\link{ictest}} }
\examples{
Aintmap(c(2,3,3,7),c(3,5,5,8),Lin=c(FALSE,TRUE,FALSE,FALSE),Rin=c(TRUE,FALSE,TRUE,FALSE))
}
\keyword{misc}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
