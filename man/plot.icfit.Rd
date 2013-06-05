\name{plot.icfit}
\alias{plot.icfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot icfit object }
\description{
Plots either the survival distributions, the cumulative distributions, or a transformation of the cumulative distributions,
from an \code{icfit} object. If there is more than one strata,
all strata will be plotted. Note that for interval censored data, the changes in the 
NPMLE of the survival function usually do not occur at unique points but occur within some interval
where any of an infinite number of curves will maximize the likelihood. We show those intervals 
were the NPMLE is indeterminate as a gray rectangle.
}
\usage{
\method{plot}{icfit}(x,XLAB="time",YLAB=NULL,COL=gray((8:1)*.1),LTY=1:9,LEGEND=NULL,
    XLEG=NULL,YLEG=NULL,shade=TRUE,dtype="survival",
    dlink=function(x){log(-log(1-x))}, xscale=1, yscale=1, conf.int=NULL,
    estpar=list( lty=NULL, lwd=1, col=gray(0)),
    cipar=list( lty=1:9, lwd=1, col=gray(0.8)),
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an icfit object, see \code{\link{icfit}}}
  \item{XLAB}{x label}
  \item{YLAB}{y label, if NULL label matches dtype}
  \item{COL}{a vector representing color of rectangles of indeterminate NPMLE, COL[i] used for ith strata}
  \item{LTY}{a vector for lty values for lines, LTY[i] used for ith strata}
  \item{LEGEND}{logical value, include legend or not, if NULL set to TRUE only if number of strata>1}
  \item{XLEG}{x location for legend, if NULL then gives maximum of 0 and minimum time from intmap}
  \item{YLEG}{y location for legend}
  \item{shade}{logical, should the rectangles of indeterminate NPMLE be colored?}
  \item{dtype}{type of distribution plotted, one of 'survival', 'cdf' or 'link' (see details)}
  \item{dlink}{link function when dtype='link' (see details)}
\item{xscale}{a numeric value used to multiply the labels on the x axis. So if the data are in days, then a value of 1/365.25 would give labels in years.}
\item{yscale}{a numeric value used to multiply the labels on the y axis. A value of 100, for instance, would be used to give a percent scale. As in the survival package, only the labels are changed, 
        not the actual plot coordinates}
\item{conf.int}{logical, should confidence intervals be plotted? NULL plots them if they are present in x object and gives no errors if they are not.}
\item{estpar}{list of par arguments for the estimated distribution lines. If lty=NULL uses LTY argument, otherwise ignores LTY (for backward compatability)}
\item{cipar}{list of par arguments for the confidence interval lines}
  \item{\dots}{other arguments passed to the plot function}
}
\details{
Turnbull (1976) noted that the NPMLE was not unique within a certain set of intervals. We represent that non-uniqueness 
using colored rectangles when \code{shade=TRUE}. The option \code{shade=TRUE} is not supported when dtype="link".

The option \code{dtype="cdf"} plots the cumulative distribution function.

When there are several strata, different types of weighted logrank-type tests (see \code{\link{ictest}}) may be derived from score statistics under the 
grouped continuous model with error distribution known. To test which test is appropriate, one may plot the cumulative distribution for 
each stratum transformed by the inverse of the proposed error distribution  (see Fay, 1996).  These are plotted with dtype="link" where dlink is the link function which 
transforms the cdf. The "wmw" scores correspond to \code{dlink=qlogis}, the "logrank2" scores correspond to the default complementary log-log  dlink, 
and the "normal" scores correspond to \code{dlink=qnorm}.
}
\value{
   Returns a list of arguments for the legend. Values are x,y, legend, fill, lty. See \code{\link{legend}} help.
}
\references{ 
    Fay, MP (1996). Rank invariant tests for interval censored data 
    under the grouped continuous model. Biometrics. 52: 811-822.

    Turnbull, B.W. (1976) The empirical distribution function 
    with arbitrarily grouped, censored and truncated data. 
    J. R. Statist.  Soc. B 38, 290-295.
 }
%\author{ ~~who you are~~ }
\note{An object of class 'icsurv' from the \code{Icens} package can use this plot function by redefining its class to 'icfit' and 'plot.icfit' will work on it.}
\seealso{\code{\link{icfit}}}
\examples{
data(bcos)
fit1<-icfit(Surv(left,right,type="interval2")~treatment,data=bcos)
summary(fit1)
plot(fit1)
 }
\keyword{hplot}
