\name{plot.icfit}
\alias{plot.icfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot icfit object }
\description{
Creates survival curves for an \code{icfit} object. If there is more than one strata,
all strata will be fit. Note that for interval censored data usually the changes in the 
NPMLE of the survival function do not occur at unique points but occur within some interval
where any of an infinite number of curves will maximize the likelihood. We show those intervals 
were the NPMLE is indeterminate as a gray rectangle.
}
\usage{
\method{plot}{icfit}(x,XLAB="time",YLAB="Survival",COL=gray((8:1)*.1),LTY=1:9,LEGEND=NULL,
    XLEG=NULL,YLEG=.1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an icfit object, see \code{\link{icfit}}}
  \item{XLAB}{x label}
  \item{YLAB}{y label}
  \item{COL}{a vector representing color of rectangles of indeterminate NPMLE, COL[i] used for ith strata}
  \item{LTY}{a vector for lty values for lines, LTY[i] used for ith strata}
  \item{LEGEND}{logical value, include legend or not, if NULL set to TRUE only if number of strata>1}
  \item{XLEG}{x location for legend, if NULL then gives maximum of 0 and minimum time from intmap}
  \item{YLEG}{y location for legend}
  \item{\dots}{other arguments passed to the plot function}
}
%\details{}
\value{
   Returns a list of arguments for the legend. Values are x,y, legend, fill, lty. See \code{\link{legend}} help.
}
%\references{ ~put references to the literature/web site here ~ }
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