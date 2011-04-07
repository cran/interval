\name{initcomputeMLE}
\alias{initcomputeMLE}
\alias{initEMICM}
\title{functions to calculate initial NPMLE of the distribution}
\description{
The function \code{\link{icfit}} calculates the NPMLE of a distribution for interval censored data using an E-M algorithm 
with polishing and checking the Kuhn-Tucker conditions (see \code{\link{icfit}} help details). 
It allows functions for the initfit option in order to calculate the starting value of the distribution in the E-M algorithm. 
Because \code{\link{icfit}} checks the Kuhn-Tucker conditions, we can try functions without doing extensive quality control,
since if the starting distribution is not close to the true NPMLE the only downside is a slower convergence. But if the 
initfit function is the true NPMLE then convergence happens on the first iteration. Functions must input 5 objects, L,R, Lin, Rin, 
and A, but need not use all of them.
}
\usage{
initcomputeMLE(L,R,Lin,Rin,A=NULL,max.inner=10,max.outer=1000,tol=1e-9)
initEMICM(L=NULL,R=NULL,Lin=NULL,Rin=NULL,A=NULL,maxiter=1000,tol=1e-7)
}
\arguments{
  \item{L}{numeric vector of left endpoints of censoring interval (equivalent to first element of Surv when type='interval2', see \code{\link{icfit}} details)}
  \item{R}{numeric vector of right endpoints of censoring interval (equivalent to second element of Surv function when type='interval2', see \code{\link{icfit}} details)}
  \item{Lin}{logical vector, should L be included in the interval? (see \code{\link{icfit}} details)}
  \item{Rin}{logical vector, should R be included in the interval? (see \code{\link{icfit}} details)} 
  \item{A}{clique matrix}
  \item{max.inner}{see \code{\link[MLEcens]{computeMLE}} }
  \item{max.outer}{see \code{\link[MLEcens]{computeMLE}} }
  \item{tol}{see either \code{\link[MLEcens]{computeMLE}} or  \code{\link[Icens]{EMICM}} }
  \item{maxiter}{see \code{\link[Icens]{EMICM}} }
}
\details{
In order to work correctly within \code{\link{icfit}} the function should output a list with at least a 'pf' element 
giving the estimated mass of the distribution for a series of intervals. 
Further, if an 'intmap' element is included (describing the series of intervals) it will be used by \code{\link{icfit}}.
}
\value{
The function \code{initcomputeMLE} outputs an icfit object with 'pf' and 'intmap' values and some other values
defined in  the help for \code{\link[MLEcens]{computeMLE}}.

The function \code{initEMICM} outputs an icsurv object with a 'pf' element but no 'intmap' element, in addition to 
some other values defined in the help for \code{\link[Icens]{EMICM}}. 

Here we define pf and intmap:
  \item{pf}{vector of estimated probabilities of the distribution}
  \item{intmap}{2 by k matrix, where the ith column defines an interval 
        corresponding to the probability, pf[i]} 
}
\section{Warning}{
In rare cases the  \code{\link[MLEcens]{computeMLE}} function (and hence the  \code{initcomputeMLE} function) can cause R to crash 
(at least for version 0.1-3 of the \code{MLEcens} package).
}
%\references{  }
\author{
The wrappers for the functions were written by M. Fay, but the real work are the calculation engines:

The calculation engine for \code{initcomputeMLE} is \code{\link[MLEcens]{computeMLE}} and was written by Marloes Maathuis, with  
part of the code for the optimization step is adapted from code that was written by Piet Groeneboom.

The calculation engine for \code{initEMICM} is \code{\link[Icens]{EMICM}} and was written by Alain Vandal and Robert Gentleman 
}
\seealso{ \code{\link{icfit}},   \code{\link[MLEcens]{computeMLE}},  \code{\link[Icens]{EMICM}} }
\examples{
## If you want speed and trust the MLEcens package, then there is no need to use icfit at all
## (but the convergence checks in icfit do not take much additional time)
data(bcos)
fit<-initcomputeMLE(bcos$left,bcos$right)
summary(fit)
plot(fit)
}
\keyword{survival}
