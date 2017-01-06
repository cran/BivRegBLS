\name{Aromatics}
\alias{Aromatics}
\docType{data}
\title{
Aromatics petroleum data
}
\description{
Aromatics measurements in light and medium petroleum by HPLC and GC MS.}
\usage{data(Aromatics)}
\format{
  A data frame with 35 observations on the following 8 variables:
  \describe{
    \item{\code{Sample}}{ a factor with the sample type}
    \item{\code{Type}}{ a factor with the following levels: HD (Heavy Diesel), LD (Light Diesel), LGO (Light Gas Oil), MGO (Medium Gas Oil).}
    \item{\code{HPLCmono}}{ a numeric vector with the monoaromatics measurements by HPLC.}
    \item{\code{GCMSmono}}{ a numeric vector with the monoaromatics measurements by GC MS.}
    \item{\code{HPLCdi}}{ a numeric vector with the diaromatics measurements by HPLC.}
    \item{\code{GCMSdi}}{ a numeric vector with the diaromatics measurements by GC MS.}
    \item{\code{HPLCtri}}{ a numeric vector with the triaromatics measurements by HPLC.}
    \item{\code{GCMStri}}{ a numeric vector with the triaromatics measurements by GC MS.}
  }
}

\source{
C-A B Ferrer, B M Celis, A B Velandia, Development of a methodology to determine the aromatic structural distribution in light and medium petroleum fractions by HPLC. Cienc. Tecnol. Futuro, 2006; 3 (2), 149-162.}
\references{
Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems, 2014; 134:123-139.}
\examples{
data(Aromatics)
str(Aromatics)
head(Aromatics)
}
\keyword{datasets}
\keyword{Aromatics petroleum data}