\name{summary.icfit}
\alias{summary.icfit}
\alias{print.icfit}
\alias{[.icfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for icfit objects}
\description{
The print method prints as a list, except the A (clique) matrix. The summary method prints the masses an associated maps for the fit. The 
[ method allows picking out of specific fits for individual elements of the factor when the right hand side of the formula in \code{icfit} 
was a factor. 
}
\usage{

\method{summary}{icfit}(object, digits=4,  \dots)

\method{print}{icfit}(x, \dots)

\method{[}{icfit}(x, i)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an icfit object}
  \item{x}{an icfit object}
  \item{digits}{number of digits for rounding results}
  \item{i}{scalar integer to pick ith strata}
  \item{\dots}{arguments to be passed}
}
%\details{}
%\value{}
%\author{ Michael P. Fay}
\seealso{ \code{\link{ictest}} }
\examples{
data(bcos)
icout<-icfit(Surv(left,right,type="interval2")~treatment, data=bcos)
print(icout)
summary(icout)
icout[1]
}
\keyword{survival}
