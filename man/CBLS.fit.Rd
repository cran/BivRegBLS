\name{CBLS.fit}
\alias{CBLS.fit}
\title{
Fit a Correlated Bivariate Least Square regression (CBLS): estimates table}
\description{
Estimate the Correlated Bivariate Least Square regression with replicated data in a (M,D) plot (Bland-Altman) where M=(X+Y)/2 and D=Y-X, provide the estimates table.}
\usage{
CBLS.fit(data = NULL, xcol = 1, ycol = 2, var.x = NULL, var.y = NULL,
     ratio.var = NULL, conf.level = 0.95)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X column(s) or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y column(s) or a character vector with the column names.}
  \item{var.x}{ a numeric variable for the variance of the measurement error of device X if known.}
  \item{var.y}{ a numeric variable for the variance of the measurement error of device Y if known.}
  \item{ratio.var}{ a numeric value for \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}, the ratio of the measurement error variances (Y over X) if known.}
  \item{conf.level}{ a numeric value for the confidence level (expressed between 0 and 1).}
}
\details{
The data argument is mandatory. If the data are unreplicated, then the measurement error variances must be given or their ratio (\ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}) in order to calculate the correlation, \ifelse{html}{\out{&rho;_{MD}}}{\eqn{\rho_{MD}}{rho_{MD}}}, between the measurement errors of the differences (on the Y-axis) and the measurement errors of the means (on the X-axis). The confidence level is used for the confidence intervals of the parameters (\ifelse{html}{\out{&rho;_{MD}}}{\eqn{\rho_{MD}}{rho_{MD}}}, \ifelse{html}{\out{&beta;}}{\eqn{\beta}{beta}} (slope), \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} (intercept)).
}
\value{
A table with the estimates of the intercept and the slope, standard error, confidence interval and pvalue (null hypothesis: slope = 0, intercept = 0).
}
\references{
Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.}
\author{Bernard G FRANCQ}

\seealso{
\code{\link{BLS}}, \code{\link{CBLS}}
}
\examples{
library(BivRegBLS)
data(SBP)
# CBLS regression on replicated data
res1=CBLS.fit(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
# CBLS regression on unreplicated data with measurement error variances previously estimated
res2=CBLS.fit(data=SBP,xcol=c("J1"),ycol="S1",var.x=80,var.y=50)
}
\keyword{CBLS}
\keyword{Correlated Bivariate Least Square}
\keyword{regression}
\keyword{confidence interval}
\keyword{confidence region}
\keyword{ellipse}
\keyword{confidence bands}
\keyword{predictive interval}
