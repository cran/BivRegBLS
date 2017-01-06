\name{DR}
\alias{DR}
\title{Deming Regression}
\description{
Estimate the Deming Regression (DR) with unreplicated or replicated data.}
\usage{
DR(data = NULL, xcol = 1, ycol = 2, ratio.var = NULL, conf.level = 0.95)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X column(s) or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y column(s) or a character vector with the column names.}
  \item{ratio.var}{ a numeric value for \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} the ratio of the measurement error variances (Y over X) if known. Otherwise, it may be estimated with replicated data.}
  \item{conf.level}{ a numeric value for the confidence level (expressed between 0 and 1).}
}
\details{
The BLS regression is more general and includes the Deming Regression. The BLS regression provides more results and should, therefore, be used instead of DR.}
\value{
A list including the following elements:
  \item{Ellipse.DR}{ a two columns matrix with the coordinates of the joint confidence interval (confidence region, ellipse) for the parameters (\ifelse{html}{\out{&beta;}}{\eqn{\beta}{beta}}, \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}}).}
  \item{Estimate.DR}{ a table (data frame) with the estimates of the intercept and the slope, standard error, confidence interval and pvalue (null hypothesis: slope = 1, intercept = 0). The exact confidence interval for the slope is also given.}
}
\references{
Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems 2014; 134:123-139.
\cr Tan CY, Iglewicz B. Measurement-methods comparisons and linear statistical relationship. Technometrics, 1999; 41(3):192-201.}
\author{Bernard G FRANCQ}
\seealso{
\code{\link{BLS}}
}
\examples{
library(BivRegBLS)
data(SBP)
res.DR=DR(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
res.DR$Estimate.DR
data(Aromatics)
res.DR=DR(data=Aromatics,xcol=3,ycol=4,ratio.var=2)
}
\keyword{Deming Regression}
\keyword{confidence interval}
\keyword{confidence region}
\keyword{ellipse}
\keyword{confidence bands}
