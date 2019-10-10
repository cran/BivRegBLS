\name{FullCIs.XY}
\alias{FullCIs.XY}
\title{
Confidence Intervals from OLSv to OLSh by DR and BLS}
\description{
Estimate the vertical and horizontal Ordinary least Square regressions, and several 'intermediate' Deming Regression (DR) and Bivariate Least Square (BLS), in a (X,Y) plot. The OLSv assumes no error on the X axis (\ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}=Infinity), while the OLSh assumes no error on the Y axis (\ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}=0). These two regressions are therefore 'extreme' regressions, while DR and BLS assume errors on both axes.}
\usage{
FullCIs.XY(data = NULL, xcol = 1, ycol = 2,
           conf.level = 0.95, npoints = 1000, nlambdas = 13)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X column(s) or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y column(s) or a character vector with the column names.}
  \item{conf.level}{ a numeric value for the confidence level (expressed between 0 and 1).}
  \item{npoints}{ an integer (at least 10) for the number of points to smooth the hyperbolic curves.}
  \item{nlambdas}{an integer for the number of intermediate DR and BLS regressions (between the OLSv and OLSh).}
}
\details{
The data argument is mandatory. This function is especially useful for unreplicated data with unknown \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} (the ratio of the measurement error variances), as it calculates all the potential solutions from OLSv to OLSh. The different estimated regression lines are provided with the different confidence intervals.}
\value{
A CIs.XY class object, a list including the following elements:
  \item{Data.means}{ a table with the X and Y data (means of the replicated data if replicated).}
  \item{Ellipses.CB}{ an array of dimension [\code{npoints}, 2 (intercept and slope), \code{nlambdas} + 2] with the coordinates of all the joint confidence intervals (confidence region, ellipses) from OLSv to OLSh.}
  \item{Slopes}{ a table (\code{nlambdas} + 2 rows) with all the slopes estimates from OLSv to OLSh including \code{nlambdas} intermediate values, the exact and approximate confidence intervals and pvalue (slope = 1).}
  \item{Intercepts}{ a table (\code{nlambdas} + 2 rows) with all the intercepts estimates from OLSv to OLSh including \code{nlambdas} intermediate values, the approximate confidence intervals and pvalue (intercept = 0).}
  \item{Joints}{ a table (\code{nlambdas} + 2 rows) with all the pvalues of the joint hypothesis (slope = 1 and intercept = 0) from OLSv to OLSh including intermediate \code{nlambdas} values.}
  \item{Hyperbolic.intervals}{ an array of dimension [npoints, 6 (X values, Y predictions, confidence interval and confidence bands),nlambdas+2] with the hyperbolic confidence intervals and confidence bands from OLSv to OLSh including \code{nlambdas} intermediate values.}
}
\references{
Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems, 2014; 134:123-139.
\cr Francq BG, Govaerts BB. Hyperbolic confidence bands of errors-in-variables regression lines applied to method comparison studies. Journal de la Societe Francaise de Statistique 2014; 155(1):23-45.}
\author{Bernard G FRANCQ}
\seealso{
\code{\link{FullCIs.MD}}, \code{\link{BLS}}
}
\examples{
library(BivRegBLS)
data(Aromatics)
res.full=FullCIs.XY(data=Aromatics,xcol=3,ycol=4)
}
\keyword{OLSv}
\keyword{OLSh}
\keyword{BLS}
\keyword{unknown ratio lambda}
