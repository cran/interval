\name{interval-internal}
\alias{icScoreTest}
\alias{icWSR}
\alias{scoresFromFit}
\alias{wsrHLYpclt}
\alias{wsrMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal functions
}
\description{
Functions mostly called by \code{\link{ictest}}. Not to be directly called by user.
}
\usage{
icScoreTest(icFIT, group, scores, alternative = "two.sided", tol.svd = 10^-8)
icWSR(fit, group, scores, alternative, type, control)
wsrHLYpclt(A, p, x, group, alternative, nwsr, type, RHO)
wsrMC(A, p, x, group, alternative, nwsr, np, digits = 12)
scoresFromFit(icFIT, scores, dqfunc=NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{icFIT}{icfit object from all the data}
  \item{group}{group vector, either numeric or character}
  \item{scores}{either 'logrank1', 'logrank2', or 'wmw'}
  \item{alternative}{either 'less','greater','two.sided', or 'two.sidedAbs'}
  \item{tol.svd}{small number, values less than that set to zero in svd}
  \item{fit}{icfit object from all the data}
   \item{type}{either 'wsr.pclt','wsr.HLY', or 'wsr.mc'}
\item{control}{list passing options to function, see \code{\link{mControl}} }
 \item{nwsr}{number of within subject resamples}
\item{np}{number of permutation replications within each wsr}
\item{digits}{for ties, treat values as tied when matched to this number of digits}
\item{A}{A matrix from icfit object on all data}
\item{p}{pf vector from icfit object on all data}
\item{x}{vector of scores associated with each value of p}
\item{dqfunc}{function, density of the results of the quantile function, e.g., function(x){dnorm(qnorm(x))} (used for 'general' scores, 
see \code{\link{ictest}})}
}
%\details{%%  ~~ If necessary, more details than the description above ~~}
%\value{}
%\references{}
%\author{}
%\note{}
%\seealso{}
%\examples{}
\keyword{internal}
