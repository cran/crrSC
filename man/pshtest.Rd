\name{psh.test}
\alias{psh.test}
\title{Goodness-of-fit test for proportional subdistribution hazards model}
\description{This Goodness-of-fit test  proposed a modified weighted Schoenfeld residuals to test the proportionality of subdistribution hazards for the Fine and Gray model
}
\usage{
psh.test(time, fstatus, z, D=c(1,1), tf=function(x) cbind(x,x^2), init)
}
\arguments{
\item{time}{vector of failure times}
\item{fstatus}{  failure status =0 if censored}
\item{z}{covariates}
 \item{D}{components of z that are tested for time-varying effect}
 \item{tf}{functions of t for z being tested on the same location}
\item{init}{initial values of regression parameters (default=all 0)}
}

\value{Returns a data.frame with percentage of cens, cause 1, Test Statistic,  d.f. ,p-value
  
}
\details{The proposed score test employs Schoenfeld residuals adapted to competing risks data. The form of the
test is established assuming that the non-proportionality arises via time-dependent coefficients in the
Fine-Gray model, similar to the test of Grambsch and Therneau.
}

\references{
Zhou B, Fine JP, Laird, G. (2013). Goodness-of-fit test for proportional subdistribution hazards mode.
Statistics in Medicine. In Press.
}


\author{Bingqing Zhou, \email{bingqing.zhou@yale.edu}}
\examples{
data(bce)
attach(bce)
lognodes <- log(nnodes)
Z1 <- cbind(lognodes, tsize/10, age, trt)  
# trt = 0 if placebo, = 0 treatment
# testing for linear time varying effect of trt
psh.test(time=time, fstatus=type, z=Z1, D=c(0,0,0,1), tf=function(x) x)
# testing for quadratic time varying effect of trt
psh.test(time=time, fstatus=type, z=Z1, D=c(0,0,0,1), tf=function(x) x^2)
# testing for log time varying effect of trt
psh.test(time=time, fstatus=type, z=Z1, D=c(0,0,0,1), 
tf=function(x) log(x))
# testing for both linear and quadratic time varying effect of trt
psh.test(time=time, fstatus=type, z=Z1, 
D=matrix(c(0,0,0,1,0,0,0,1), 4,2),  tf=function(x) cbind(x,x^2))

}

