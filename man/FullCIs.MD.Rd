\name{FullCIs.MD}
\alias{FullCIs.MD}
\title{
Confidence Intervals by CBLS with all potential solutions}
\description{
Estimate the Correlated-Bivariate-Least Square regression (CBLS) for all potential solutions from \ifelse{html}{\out{&rho;MD}}{\eqn{\rho_{MD}}{rhoMD}}=-1 to \ifelse{html}{\out{&rho;MD}}{\eqn{\rho_{MD}}{rhoMD}}=1, in a (M,D) plot. This function is analogous to \code{\link{FullCIs.XY}} estimating all the BLS regressions from OLSv to OLSh in a (X,Y) plot.}
\usage{
FullCIs.MD(data = NULL, xcol = 1, ycol = 2,
           conf.level = 0.95, npoints = 1000, nlambdas = 13)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X column(s) or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y column(s) or a character vector with the column names.}
  \item{conf.level}{ a numeric value for the confidence level (expressed between 0 and 1).}
  \item{npoints}{ a numeric variable for the number of points to smooth the hyperbolic curves.}
  \item{nlambdas}{a numeric variable for the number of intermediate CBLS regressions (the two extreme CBLS are calculated by default).}
}
\details{
The data argument is mandatory. This function is especially useful for unreplicated data with unknown \ifelse{html}{\out{&rho;MD}}{\eqn{\rho_{MD}}{rhoMD}} (related to \ifelse{html}{\out{&lambda;XY}}{\eqn{\lambda_{XY}}{lambdaXY}}, the ratio of the measurement error variances), as it calculates all the potential solutions. The different estimated regression lines are provided with the different confidence intervals.}
\value{
A CIs.MD class object, a list including the following elements:
  \item{Data.means}{ a table with the X and Y data (means of the replicated data if replicated), their means and differences.}
  \item{Ellipses.CB}{ an array of dimension [\code{npoints}, 2 (intercept and slope), \code{nlambdas} + 2] with the coordinates of all the joint confidence intervals (confidence region, ellipses) from \ifelse{html}{\out{&rho;MD}}{\eqn{\rho_{MD}}{rhoMD}} = -1 to \ifelse{html}{\out{&rho;MD}}{\eqn{\rho_{MD}}{rhoMD}} = 1.}
  \item{Slopes}{ a table (\code{nlambdas} + 2 rows) with all the slopes estimates and their approximate confidence intervals and pvalue (slope = 0).}
  \item{Intercepts}{ a table (\code{nlambdas} + 2 rows) with all the intercepts estimates and their approximate confidence intervals and pvalue (intercept = 0).}
  \item{Joints}{ a table (\code{nlambdas} + 2 rows) with all the pvalues of the joint hypothesis (slope = 0 and intercept = 0).}
  \item{Hyperbolic.intervals}{ an array of dimension [\code{npoints}, 6 (X values, Y predictions, confidence interval and confidence bands), \code{nlambdas} + 2] with the hyperbolic confidence interval and confidence bands.}
}
\references{
Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.
}
\author{
Bernard G FRANCQ}
\seealso{
\code{\link{FullCIs.XY}}, \code{\link{CBLS}}
}
\examples{
library(BivRegBLS)
data(Aromatics)
res.full=FullCIs.MD(data=Aromatics,xcol=3,ycol=4)
}
\keyword{CBLS}
\keyword{unknown ratio lambda}
