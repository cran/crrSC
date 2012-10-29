\name{cdata}
\Rdversion{1.1}
\alias{cdata}
\docType{data}
\title{Clustered competing risks simulated  data}
\description{sample of 200 observations}

\usage{data(cdata)}
\format{
A data frame with 200 observations and the following 4 variables. 
Simulation is detailed on the paper Competing Risk Regression for clustered data. 
Zhou, Fine, Latouche, Labopin. 2011. In Press. Biostatistics.
  \describe{

    \item{\code{ID}}{Id  of cluster, each cluster is of size 2}
    \item{\code{ftime}}{Event time}
    \item{\code{fstatus}}{Event type. 0=censored, 1 , 2}
    \item{\code{z}}{a binary covariate with P(z=1)=0.5}
    
   
  }
}

\examples{data(cdata)}
\keyword{datasets}
\keyword{crrSC}
