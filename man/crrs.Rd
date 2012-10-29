\name{crrs}
\alias{crrs}
\title{
Competing Risks Regression for Stratified Data
}
\description{Regression modeling of subdistribution hazards for stratified right censored data 

Two types of stratification are addressed :
 Regularly stratified: small number of large groups (strata) of subjects
 Highly stratified: large number of small groups (strata) of subjects

}
\usage{

crrs(ftime, fstatus, cov1, cov2, strata, 
tf, failcode=1, cencode=0,
ctype=1, 
subsets, na.action=na.omit, 
gtol=1e-6, maxiter=10,init)
}
\arguments{
\item{strata}{stratification covariate
}
\item{ctype}{1 if estimating censoring dist within strata (regular stratification), 
2 if estimating censoring dist across strata (highly stratification)
}

\item{ftime}{
vector of failure/censoring times
}
\item{fstatus}{
vector with a unique code for each failure type and a separate code for
censored observations 
}
\item{cov1}{
matrix (nobs x ncovs) of fixed covariates (either cov1, cov2, or both
are required)
}
\item{cov2}{
matrix of covariates that will be multiplied by functions of time; 
if used, often these covariates would also appear in cov1
to give a prop hazards effect plus a time interaction
}
\item{tf}{
functions of time.  A function that takes a vector of times as
an argument and returns a matrix whose jth column is the value of 
the time function corresponding to the jth column of cov2 evaluated
at the input time vector.  
At time \code{tk}, the model includes the term \code{cov2[,j]*tf(tk)[,j]} as a covariate.
}

\item{failcode}{
code of fstatus that denotes the failure type of interest
}
\item{cencode}{
code of fstatus that denotes censored observations
}
\item{subsets}{
  a logical vector specifying a subset of cases to include in the
  analysis
}
\item{na.action}{
  a function specifying the action to take for any cases missing any of
  ftime, fstatus, cov1, cov2, cengroup, or subset.
}
\item{gtol}{
  iteration stops when a function of the gradient is \code{< gtol}
}
\item{maxiter}{
maximum number of iterations in Newton algorithm (0 computes
scores and var at \code{init}, but performs no iterations)
}
\item{init}{
  initial values of regression parameters (default=all 0)
}
}
\value{Returns a list of class crr, with components (see crr for details)
  \item{$coef}{the estimated regression coefficients}
  \item{$loglik}{log pseudo-liklihood evaluated at \code{coef}}
\item{$score}{derivitives of the log pseudo-likelihood evaluated at \code{coef}}
\item{$inf}{-second derivatives of the log pseudo-likelihood}
\item{$var}{estimated variance covariance matrix of coef}
\item{$res}{matrix of residuals}
\item{$uftime}{vector of unique failure times}
\item{$bfitj}{jumps in the Breslow-type estimate of the underlying
  sub-distribution cumulative hazard (used by predict.crr())}
\item{$tfs}{the tfs matrix (output of tf(), if used)}
\item{$converged}{TRUE if the iterative algorithm converged}
\item{$call}{The call to crr}
\item{$n}{The number of observations used in fitting the model}
\item{$n.missing}{The number of observations removed from the input data
  due to missing values}
\item{$loglik.null}{The value of the log pseudo-likelihood when all the
  coefficients are 0}
}
\details{Fits the stratified extension of the Fine and Gray model (2011).  
This model directly assesses the effect of covariates on the subdistribution 
of a particular type of failure in a competing risks setting.  
}

\references{
Zhou B, Latouche A, Rocha V, Fine J. (2011)  Competing Risks Regression for Stratified Data.
Biometrics, 67(2):661-70.
}

\seealso{cmprsk}
\author{Bingqing Zhou, \email{bingqing.zhou@yale.edu}}
\examples{
##
#using fine and gray model
#crr(ftime=center$ftime, fstatus=center$fstatus, 
#cov1=cbind(center$fm,center$cells))
#
# High Stratification: ctype=2
# Random sub-sample 
data(center)
cov.test<-cbind(center$fm,center$cells)
crrs(ftime=center[,1],fstatus=center[,2],
cov1=cov.test,
strata=center$id,ctype=2)

}
\keyword{survival}
