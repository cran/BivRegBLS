\name{SBP}
\alias{SBP}
\docType{data}
\title{
Systolic blood pressure data
}
\description{
Systolic blood pressure measured by two devices on 85 patients: 3 times with a manual device by 2 operators and 3 times with a semi automatic device.
}
\usage{data(SBP)}
\format{
  A data frame with 85 observations and 10 variables.
  \describe{
    \item{\code{Subject}}{a numeric vector for the patient IDs.}
    \item{\code{J1}}{a numeric vector: the 1st measures obtained by operator J.}
    \item{\code{J2}}{a numeric vector: the 2nd measures obtained by operator J.}
    \item{\code{J3}}{a numeric vector: the 3rd measures obtained by operator J.}
    \item{\code{R1}}{a numeric vector: the 1st measures obtained by operator R.}
    \item{\code{R2}}{a numeric vector: the 2nd measures obtained by operator R.}
    \item{\code{R3}}{a numeric vector: the 3rd measures obtained by operator R.}
    \item{\code{S1}}{a numeric vector: the 1st measures obtained by the semi automatic device.}
    \item{\code{S2}}{a numeric vector: the 2nd measures obtained by the semi automatic device.}
    \item{\code{S3}}{a numeric vector: the 3rd measures obtained by the semi automatic device.}
  }
}
\source{
Bland JM, Altman DG. Measuring agreement in method comparison studies. Statistical Methods in Medical Research, 1999; 8:135-160.}
\references{
Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems, 2014; 134:123-139.
\cr Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.
}
\examples{
data(SBP)
str(SBP)
head(SBP)
}
\keyword{datasets}
\keyword{systolic blood pressure data}
