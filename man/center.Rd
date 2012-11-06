\name{center}
\Rdversion{1.1}
\alias{center}
\docType{data}
\title{Multicenter Bone  Marrow transplantation data}
\description{Random sub sample of 400 patients}

\usage{data(center)}
\format{
  A data frame with 400 observations and the following 5 variables.
  \describe{
    \item{\code{id}}{Id of transplantation center}
    \item{\code{ftime}}{Event time}
    \item{\code{fstatus}}{Event type. 0=censored, 1=Acute or Chronic GvHD , 2=Death free of GvHD}
    \item{\code{cells}}{source of stem cells: peripheral blood vs bone marrow}
    \item{\code{fm}}{female donor to male recipient match}
   
  }
}

\examples{data(center)}
\keyword{datasets}
\keyword{crrSC}
