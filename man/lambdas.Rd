\name{lambdas}
\alias{lambdas}
\title{
Measurement error variances ratio}
\description{
Calculate the measurement error variances ratio of two devices (Y over X): \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} and \ifelse{html}{\out{&lambda;_{XY}}}{\eqn{\lambda_{XY}}{lambda_{XY}}.}
}
\usage{
lambdas(data = NULL, xcol = NULL, ycol = NULL, conf.level = 0.95)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X columns or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y columns or a character vector with the column names.}
  \item{conf.level}{ a numeric value for the confidence level.}
}
\details{
The data must be replicated to estimate the measurement error variances. If the number of replicates in X is equal to the number of replicates in Y, then \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} and \ifelse{html}{\out{&lambda;_{XY}}}{\eqn{\lambda_{XY}}{lambda_{XY}}} are equal: \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} is the ratio (Y over X) of the measurement error variances, while \ifelse{html}{\out{&lambda;_{XY}}}{\eqn{\lambda_{XY}}{lambda_{XY}}} is similar but takes also into account the number of replicates per device (nx and ny). Unbiased estimators (which is not the ratio of the two variances) for \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} and \ifelse{html}{\out{&lambda;_{XY}}}{\eqn{\lambda_{XY}}{lambda_{XY}}} are also given.}
\value{
A lambdas class object, a table with 2 rows (\ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} and \ifelse{html}{\out{&lambda;_{XY}}}{\eqn{\lambda_{XY}}{lambda_{XY}}}) and their confidence intervals and pvalues in columns (the null hypothesized value is 1).
}
\references{
Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.
\cr Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems, 2014; 134:123-139.}
\author{Bernard G FRANCQ}
\seealso{
\code{\link{desc.stat}}
}
\examples{
library(BivRegBLS)
data(SBP)
lambdas(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
}
\keyword{variances}
\keyword{ratio}
\keyword{lambda}
