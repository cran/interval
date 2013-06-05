\name{mControl}
\alias{mControl}
\title{Auxiliary for feeding parameters to different methods}
\description{
A function to create a list of arguments for \code{\link{ictest}}.
}
\usage{
mControl(cm=NULL,nmc=10^3-1,seed=1234321,digits=12,p.conf.level=.99,
    setSEED=TRUE,tol.svd=10^-8,nwsr=10^3-1,np=10^3-1,tsmethod="central")
}
\arguments{
\item{cm}{a choose(n,m) by n matrix, used if method='exact.ce', ignored otherwise}
\item{nmc}{number of Monte Carlo replications, used if method='exact.mc', ignored otherwise}
\item{seed}{value used in \code{set.seed} if method='exact.mc', or any of three wsr methods, ignored otherwise} 
\item{setSEED}{logical, set to FALSE when performing simulations}
\item{p.conf.level}{confidence level for p value estimate, used if method='exact.mc', ignored otherwise }
\item{digits}{number of digits to use in \code{\link{signif}} for precision of test statistics}
\item{tol.svd}{tolerance for use in calculating g-inverse, values less than tol.svd are set to zero, used when method='scoretest'}
\item{nwsr}{number of within subject resamples, used when method='wsr.mc', 'wsr.HLY', or 'wsr.pclt'}
\item{np}{number of permutation replications within each wsr, used when method='wsr.mc'}
\item{tsmethod}{two-sided method for exact permutation tests, either 'central' or 'abs' (see \code{\link[perm]{permControl}})}
}
\details{

When \code{cm}=NULL the resulting matrix is created by \code{\link[perm]{chooseMatrix}}, it may be optionally provided here 
only so that \code{chooseMatrix} does not need to be repeatedly called in simulations.  
Also when doing simulations (with
method='exact.mc' or any of the wsr methods), use setSEED=FALSE so that the seed 
is not reset to the same value each time you call the function.

See \code{\link[perm]{calcPvalsMC}} for description of how p.conf.level is used. 
}
\value{
An list with the arguments as components.
}
%\author{Michael Fay}
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
