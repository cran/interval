\name{icfit}
\alias{icfit}
\alias{icfit.default}
\alias{icfit.formula}
%\alias{[.icfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calculate non-parametric MLE for interval censored survival function }
\description{
This function calculates the the non-parametric maximum likelihood estimate for the distribution from interval 
censored data using the self-consistent estimator, so the associated survival distribution generalizes 
the Kaplan-Meier estimate to interval censored data. Formulas using \code{Surv} are allowed similar to \code{survfit}.
}
\usage{

\method{icfit}{formula}(formula, data, \dots)

\method{icfit}{default}(L, R,initfit =NULL, control=icfitControl(), Lin=NULL, Rin=NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{numeric vector of left endpoints of censoring interval (equivalent to first element of \code{Surv} when type='interval2', see details)}
  \item{R}{numeric vector of right endpoints of censoring interval (equivalent to second element of \code{Surv} function when type='interval2', see details)}
  \item{initfit}{an initial estimate as an object of class \code{icfit} or \code{icsurv}, or a character vector of the name of the function used to calculate the initial estimate (see details)}
  \item{control}{list of arguments for controling algorithm (see \code{\link{icfitControl}})}
  \item{Lin}{logical vector, should L be included in the interval? (see details)}
  \item{Rin}{logical vector, should R be included in the interval? (see details)} 
  \item{formula}{a formula with response a numeric vector (which assumes no censoring) or \code{Surv} object
          the right side of the formula may be 1 or a factor (which produces separate fits for each level). }
  \item{data}{an optional matrix or data frame containing the variables in the formula. By default the variables are taken from environment(formula).
}
  \item{\dots}{values passed to other functions}
}
\details{
The \code{icfit} function fits the nonparametric maximum likelihood estimate (NPMLE) of the 
distribution function for interval censored data. In the default case (when Lin=Rin=NULL) 
we assume there are n (n=length(L)) failure times, and the ith one is in the interval
between L[i] and  R[i]. The default is not to include L[i] in the interval unless L[i]=R[i],
and to include R[i] in the interval unless R[i]=Inf. When Lin and Rin are not NULL they describe 
whether to include L and R in the associated interval. If either Lin or Rin is length 1 then it is 
repeated n times, otherwise they should be logicals of length n. 

   The algorithm is basically an EM-algorithm applied to 
    interval censored data (see Turnbull, 1976); however 
    first we can define a set of intervals (called the Turnbull intervals)
    which are the only intervals where the NPMLE may change. The Turnbull intervals are also called the 
    innermost intervals, and are the result of the primary reduction (see Aragon and 
    Eberly, 1992). The starting distribution for the E-M algorithm is given by \code{initfit}, which may be either 
 (1) NULL, in which case a very simple and quick starting distribution is used (see code), (2) a character vector 
describing a function with inputs, L,R, Lin, Rin, and A, see for example \code{\link{initcomputeMLE}}, (3)
a list giving \code{pf} and \code{intmap} values, e.g., an \code{icfit} object. If option (2) is tried and results in an error then 
the starting distribution reverts to the one used with option (1).
Convergence is defined when the maximum 
    reduced gradient is less than epsilon (see \code{\link{icfitControl}}), and the 
    Kuhn-Tucker conditions are approximately met, 
    otherwise a warning will result.  (see Gentleman and 
    Geyer, 1994). There are other faster algorithms (for example see 
\code{\link[Icens]{EMICM}} in the package 
   \code{Icens}.

The output is of class \code{icfit} which is identical to the \code{icsurv} class of the 
\code{Icens} package when there is only one group for which a distribution is needed. 
Following that class, there is an \code{intmap} element which gives the bounds 
about which each drop in the NPMLE survival function can occur.  


Since the classes \code{icfit} and \code{icsurv} are so closely related, one can directly 
use of initial (and faster) fits from the \code{Icens} package as input in 
\code{initfit}. Note that when using a non-null \code{initfit}, the \code{Lin} and \code{Rin} values of the 
initial fit are ignored. Alternatively, one may give the name of the function used to calculate the initial fit. 
The function is assumed to input the transpose of the A matrix (called A in the Icens package). Options can be passed 
to initfit function as a list using the initfitOpts variable in \code{\link{icfitControl}}.

The advantage of the \code{icfit} function over those in \code{Icens} package is that it allows a call similar 
to that used in \code{survfit} of the \code{survival} package so that different groups may be 
plotted at the same time with similar calls.  

An \code{icfit} object prints as a list (see value below). A \code{print} function prints output as a list 
except suppresses printing of A matrix. A \code{summary} function prints the 
distribution (i.e., probabilities and the intervals where those 
probability masses are known to reside) for each group in the icfit object. There is also 
a plot method, see \code{\link{plot.icfit}}.

For additional references and background see Fay and Shaw (2010).
}
\value{
An object of class \code{icfit} (same as icsurv class, see details). 

There are 4 methods for this class: \code{\link{plot.icfit}}, \code{\link{print.icfit}}, \code{\link{summary.icfit}}, and \code{\link{[.icfit}}. The last method 
pulls out individual fits when the right side of the formula of the \code{icfit} call was a factor.  

A list with elements:
  \item{A}{ this is the n by k matrix of indicator functions, NULL if more than one strata, not printed by default}
  \item{strata}{a named numeric vector of numbers of observations in each strata, if one strata observation named NPMLE} 
  \item{error}{this is max(d + u - n), see Gentleman and Geyer, 1994}
  \item{numit}{number of iterations}
  \item{pf}{vector of estimated probabilities of the distribution}
  \item{intmap}{2 by k matrix, where the ith column defines an interval 
        corresponding to the probability, pf[i]}
  \item{converge}{a logical, TRUE if normal convergence}
  \item{message}{character text message on about convergence}
  \item{anypzero}{logical denoting whether any of the Turnbull intervals were set to zero}
}
\references{    Aragon, J and Eberly, D (1992). On convergence of convex 
    minorant algorithms for distribution estimation with 
    interval-censored data.  J. of Computational and Graphical 
    Statistics. 1: 129-140.

    Fay, MP and Shaw, PA (2010). Exact and Asymptotic Weighted Logrank Tests for Interval Censored Data: The interval R package.
    Journal of Statistical Software. \url{http://www.jstatsoft.org/v36/i02/}. 36 (2):1-34.

    Gentleman, R. and Geyer, C.J. (1994). Maximum likelihood 
    for interval censored data:consistency and computation. 
    Biometrika, 81, 618-623.

    Turnbull, B.W. (1976) The empirical distribution function 
    with arbitrarily grouped, censored and truncated data. 
    J. R. Statist.  Soc. B 38, 290-295.
}
\author{ Michael P. Fay}
\seealso{ \code{\link{ictest}},  \code{\link[Icens]{EMICM}} }
\examples{
data(bcos)
icout<-icfit(Surv(left,right,type="interval2")~treatment, data=bcos)
plot(icout)
## can pick out just one group
plot(icout[1])
}
\keyword{survival}
