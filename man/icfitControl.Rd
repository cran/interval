\name{icfitControl}
\alias{icfitControl}
\title{Auxiliary for controling icfit }
\description{
A function to create a list of arguments for \code{\link{icfit}}.
}
\usage{
icfitControl(epsilon = 1e-06, maxit = 10000, initfitOpts=NULL, conf.level=.95, B=200, 
    confMethod="modboot",seed=19439101, timeEpsilon=1e-06, timeMessage=TRUE) 
}
\arguments{ 
\item{epsilon}{The minimum error for convergence purposes. The 
	   EM algorithm stops when error \eqn{<}{<} epsilon, where 
	   error is the maximum of the reduced gradients 
           (see Gentleman and Geyer, 1994)}
  \item{maxit}{maximum number of iterations of the EM algorithm}
    \item{initfitOpts}{named list of options for initfit function if initfit is function name}
\item{conf.level}{level for confidence interval, used if conf.int=TRUE}
\item{B}{number of bootstrap replications for conf.int=TRUE, must be at least 11}
\item{confMethod}{method for confidence intervals, must be "modboot"}
\item{seed}{random seed for bootstrap, if NULL no call to set.seed}
\item{timeEpsilon}{small number for adding or subtracting from time for drawing confidence interval lines}
\item{timeMessage}{logical, print estimate of how long modified bootstrap confidence intervals will take to calculate?}  
}
\details{
There is only one option for the confMethod now. The confMethod argument is only needed for future versions if there is another confidence interval method option.
For a description of the modified bootstrap method see \code{\link{icfit}}. 
}
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
