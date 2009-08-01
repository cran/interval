\name{ictest}
\alias{ictest}
\alias{wlr_trafo}
\alias{ictest.default}
\alias{ictest.formula}
\alias{wlr_trafo.default}
\alias{wlr_trafo.Surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{do logrank or Wilcoxon type tests on interval censored data}
\description{
 The \code{ictest} function performs several different tests for 
interval censored data, and the \code{wlr_trafo} function gives the 
associated transformation to scores. The default tests use the permutation 
form of the test (the method option tells whether it is asymptotic or exact), but the score test 
form and multiple imputation form are supported. 
The 3 different scores give three different tests 
that generalize to interval censored data either the 
Wilcoxon-Mann-Whitney test
(scores="wmw") or the logrank test (scores="logrank1" or 
scores="logrank2") (see details). 


The function calls the \code{\link{icfit}} function, if an icfit object is not provided.
}
\usage{

\method{ictest}{default}(L, R, group,  
    scores = c("logrank1","logrank2","wmw"),
    rho=NULL,
    alternative= c("two.sided", "less", "greater","two.sidedAbs"),   
    icFIT=NULL,
    initfit=NULL, 
    icontrol=icfitControl(),
    exact=NULL,
    method=NULL,
    methodRule=methodRuleIC1,
    mcontrol=mControl(),
    Lin=NULL,
    Rin=NULL,\dots)

\method{ictest}{formula}(formula, data, subset, na.action, \dots)




\method{wlr_trafo}{default}(x, R=NULL, 
    scores = c("logrank1", "logrank2", "wmw"), 
    icFIT = NULL, initfit = NULL, control=icfitControl(),
    Lin=NULL,Rin=NULL,\dots)

\method{wlr_trafo}{Surv}(x,\dots)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{numeric vector of left endpoints of censoring interval, if R is NULL
then represents exact failure time}
  \item{R}{numeric vector of right endpoints of censoring interval}
 \item{x}{response, either a Surv object or a numeric vector representing
the left endpoint. if latter and R is NULL then x is treated as exact} 
  \item{group}{ 	a vector denoting the group for which the test 
	is desired. If group is a factor or character 
	then a k-sample test is performed, where k is 
	the number of unique values of group. If group
	is numeric then a "correlation" type test is 
	performed. If there are only two groups, both 
	methods give the same results.  }
  \item{scores}{character vector defining the scores: "logrank1" (default), "logrank2", or "wmw" (see details)}
  \item{rho}{either 0 (gives scores="logrank1"), or 1 (gives scores="wmw") (see Note)}
\item{alternative}{character giving alternative for two-sample and trend tests, K-sample should be two.sided} 
 \item{icFIT}{a precalculated icfit object for increased computation speed. This should be the icfit from 
the pooled data. Normally initfit should be used instead (see Warning)}
  \item{initfit}{an object of class icfit or icsurv, used for the initial estimate (see Warning). 
Ignored if \code{icFIT} is not null}
   \item{icontrol}{list of arguments for controling NPMLE algorithm in call to icfit (default \code{\link{icfitControl}})}
  \item{formula}{a formula with response a numeric vector (which assumes no censoring) or Surv object,
          the right side of the formula is the group variable. No strata() is allowed}
\item{data}{data frame for variables in formula}
\item{subset}{  an optional vector specifying a subset of observations to be used }
\item{na.action}{ a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action")}
\item{Surv}{ a Surv object, see \code{\link[survival]{Surv}} }
\item{exact}{a logical value, TRUE denotes exact test, ignored if method is not NULL}
\item{method}{a character value, one of 'pclt','exact.network','exact.ce','exact.mc', 'scoretest', 'wsr.HLY', 'wsr.pclt', 'wsr.mc'. If NULL method chosen by methodRule.}
\item{methodRule}{a function used to choose the method, default \code{\link{methodRuleIC1}}. (see details in \code{\link[perm]{perm}})}
 \item{mcontrol}{list of arguments for controling algorithms of different methods (see \code{\link{mControl}})}
  \item{Lin}{logical vector, should L be included in the interval? (see details)}
  \item{Rin}{logical vector, should R be included in the interval? (see details)} 
   \item{control}{list of arguments for controling NPMLE algorithm in call to icfit (default \code{\link{icfitControl}})}
  \item{\dots}{values passed to other functions}
}
\details{

The censoring in the default case (when Lin=Rin=NULL) 
assumes there are n (n=length(L)) failure times, and the ith one is in the interval
between L[i] and  R[i]. The default is not to include L[i] in the interval unless L[i]=R[i],
and to include R[i] in the interval unless R[i]=Inf. When Lin and Rin are not NULL they describe 
 whether to include L and R in the associated interval. If either Lin or Rin is length 1 then it is 
repeated n times, otherwise they should be logicals of length n. 


   The 3 different types of scores are compared in depth in Fay (1999). When 
   scores='logrank1' this gives the most commonly used logrank scores for right censored data,
   and reduces to the scores of Sun (1996) for interval censored data.  When 
   scores='logrank2' this gives the scores associated with the grouped 
   proportional hazards  model of Finkelstein (1986). When scores='wmw' this 
   gives the generalized Wilcoxon-Mann-Whitney scores. 

   For censored data 2 common likelihoods are the 
   marginal likelihood of the ranks and the 
   likelihood with nuisance parameters for the 
   baseline survival. Here we use the latter 
   likelihood (as in Finkelstein, 1986, Fay, 1996,
   and Sun, 1996). 

Because of theoretical difficulties (discussed below), the default method is to
   perform a permutation test on the scores. There are several ways to perform the permutation 
test, and the function \code{\link{methodRuleIC1}} chooses which of these ways will be used. The choice
is basically between using a permutational central limit theorem (method="pclt") or using an exact method. 
There are several algorithms for the exact method (see \code{\link[perm]{perm}} ). 
   
Another method is to perform a standard score test (method="scoretest").
It is difficult to prove the asymptotic validity of the standard
   score tests for this likelihood because the 
   number of nuisance parameters typically grows with
   the sample size and often many of 
   the parameters are equal at the nonparametric 
   MLE, i.e., they are on the boundary of the 
   parameter space (Fay, 1996). Specifically, when the score test is performed then 
   an adjustment is made so that the nuisance parameters are defined based on the data and do not approach 
the boundary of the parameter space (see Fay, 1996).  Theoretically, the score test should perform well when there are 
many individuals but few observation times, and its advantage in this situation is that it retains validity even when 
the censoring mechanism
may depend on the treatment.   

Another method is to use multiple imputation, or within subject resampling (method="wsr.HLY") (Huang, Lee, and Yu, 2008). 
This method samples interval censored observations from the nonparametric distribution, then performs the usual 
Martingale-based variance. A different possibility is to use a permutational central limit theorem variance for each
wsr (method="wsr.pclt") or use Monte Carlo replications to get an possibly exact method from each
within subject resampling (method="wsr.mc").  

   Note that when icfit and ictest are used on right 
   censored data, because of the method of estimating
   variance is different, even Sun's method does not
   produce exactly the standard logrank test results. 

   There are some typos in Appendix II of Fay (1999), 
   see the correction.

\code{\link[perm]{permControl}} 

}
\value{
	The function \code{wlr_trafo} returns only the numeric vector of scores, while 
\code{ictest} returns an object of class `ictest', which is a list with the following values. 
	
\item{scores}{
  This is a vector the same length as L and R, containing the scores
  used in the permutation test. }

\item{U}{
  The efficient score vector. When group is a factor 
  or character vector then each element of U has the 
  interpretation as the weighted sum of "observed" minus
  "expected" deaths for the group element defined by the 
  label of U. Thus negative values indicate better than average
  survival (see Fay, 1999). }

\item{N}{ number of observations in each group }

\item{method}{ full description of the test }
\item{data.name}{ description of data variables} 

\item{algorithm}{ algorithm used in permutation calculation. For example, 
`pclt', 'exact.mc', 'exact.network'}

\item{statistic}{ either the chi-square or Z statistic, or NULL for exact methods}
\item{parameter}{degrees of freedom for chi-square statistic}
\item{alternative}{alternative hypothesis}
\item{alt.phrase}{ phrase used to describe the alternative hypothesis }

\item{p.value}{p value associated with alternative}
\item{p.values}{ vector of p-values under different alternatives} 
\item{p.conf.int}{ confidence interval on p.value, for algorithm='exact.mc' only }
\item{nmc}{number of Monte Carlo replications, for algorithm='exact.mc' only }
\item{nwsr}{number of within subject resamplings, for WSR methods only}
\item{np}{number of permutation replications within each WSR, for method='wsr.mc' only}
}
\references{    Fay, MP (1996). "Rank invariant tests for interval 
   censored data under the grouped continuous model".
   Biometrics, 52: 811-822.

   Fay, MP (1999). "Comparing Several Score Tests for 
   Interval Censored Data." Statistics in Medicine,
   18: 273-285 (Correction: 1999, 18: 2681).

   Finkelstein, DM (1986). "A proportional hazards 
   model for interval censored failure time data"
   Biometrics, 42: 845-854.

   Huang, J, Lee, C, Yu, Q (2008). "A generalized log-rank test for interval-censored failure time 
   data via multiple imputation" Statistics in Medicine, 27: 3217-3226.

   Sun, J (1996). "A non-parametric test for interval
   censored failure time data with applications to 
   AIDS studies". Statistics in Medicine, 15: 1387-1395. 

}
\author{ Michael P. Fay }
\note{The \code{rho} argument gives the scores which match the scores from the 
\code{\link[survival]{survdiff}} function, so that when rho=0 then scores="logrank1",
and when rho=1 then scores="wmw". These scores will exactly match those used in survdiff,
but the function \code{survdiff} uses an asymptotic method 
based on the score test to calculate p-values, while \code{ictest} uses 
permutation methods to calculate 
the p-values, so that the p-values will not match exactly.
}
\section{Warning }{Because the input of \code{icfit} is only for saving computational time, 
no checks are made to determine if the \code{icfit} is in fact the correct one. Thus you will get 
wrong answers with no warnings if you input the wrong \code{icfit} object. The safer way to save 
computational time is to input the precalculated \code{icfit} object into \code{initfit}. When this
is done, you will get either the correct answer or a warning even when you input a bad guess for the 
initfit.} 
\seealso{ \code{\link{icfit}} }
\examples{
   ## perform a logrank-type test using the permutation form of the test
data(bcos)
testresult<-ictest(Surv(left,right,type="interval2")~treatment, scores="logrank1",data=bcos)
testresult
## perform a Wilcoxon rank sum-type test
## using asymptotic permutation variance
left<-bcos$left
right<-bcos$right
trt<-bcos$treatment
## save time by using previous fit
ictest(left,right,trt, initfit=testresult$fit, method="pclt",scores="wmw")
}
\keyword{htest }% at least one, from doc/KEYWORDS
\keyword{survival }% __ONLY ONE__ keyword per line
\keyword{nonparametric}