\name{bce}
\Rdversion{1.1}
\alias{bce}
\docType{data}
\title{Breast Cancer Data}
\description{Data Randomly Generated According To El178 clinical trial}

\usage{data(bce)}
\format{
  A data frame with 200 observations and the following 6 variables.
  \describe{
    \item{\code{trt}}{Treatment: 0=Placebo, 1=Tamoxifen}
    \item{\code{time}}{Event time}
    \item{\code{type}}{Event type. 0=censored, 1=Breast Cancer recurrence , 2=Death without recurrence }
    \item{\code{nnodes}}{Number of positive nodes}
    \item{\code{tsize}}{Tumor size}
    \item{\code{age}}{Age}   
  }
}

\examples{data(bce)}
\keyword{datasets}
\keyword{crrSC}



