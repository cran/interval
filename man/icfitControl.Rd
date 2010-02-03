\name{icfitControl}
\alias{icfitControl}
\title{Auxiliary for controling icfit }
\description{
A function to create a list of arguments for \code{\link{icfit}}.
}
\usage{
icfitControl(epsilon = 1e-06, maxit = 10000, initfitOpts=NULL) 
}
\arguments{ 
\item{epsilon}{The minimum error for convergence purposes. The 
	   EM algorithm stops when error \eqn{<}{<} epsilon, where 
	   error is the maximum of the reduced gradients 
           (see Gentleman and Geyer, 1994)}
  \item{maxit}{maximum number of iterations of the EM algorithm}
    \item{initfitOpts}{named list of options for initfit function if initfit is function name}
}
%\details{}
\references{  
   Gentleman, R. and Geyer, C.J. (1994). Maximum likelihood 
    for interval censored data:consistency and computation. 
    Biometrika, 81, 618-623.

  }

\value{
An list with the arguments as components.
}
%\author{Michael Fay}
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
