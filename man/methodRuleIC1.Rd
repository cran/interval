\name{methodRuleIC1}
\alias{methodRuleIC1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rule for determining method for ictest}
\description{
This is the default function which determines which permutation method (e.g., `pclt' or `exact.network') to use in \code{\link{ictest}}.  
}
\usage{
methodRuleIC1(x, group, exact, Nbound = c(20))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of response scores  }
  \item{group}{group membership vector}
  \item{exact}{logical, TRUE=exact method chosen, FALSE=pclt}
  \item{Nbound}{bound, if n>Nbound then method='pclt' otherwise either 'exact.mc' (for k-sample or trend) 
       or 'exact.network' (for two-sample)}
}
\details{
This function determines which of several methods will be used in 
\code{\link{ictest}}, see \code{\link[perm]{permTS}} for description of methods. 

When exact=FALSE then returns 'pclt'. When exact=TRUE then returns either 'exact.network' if the 
length(cc)<=Nbound and it is a two-sample test or 'exact.mc' otherwise. When exact=NULL and the 
length(cc)<=Nbound, then returns either 'exact.network' (for two-sample) or 'exact.mc' (for k-sample 
and trend). When exact=NULL and length(cc)>Nbound returns 'pclt'.
}
\value{
a character vector with one of the following values: "pclt","exact.network","exact.mc"
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
\seealso{ \code{\link{ictest}}}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
