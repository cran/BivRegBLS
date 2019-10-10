\name{GraphFullCIs.XY}
\alias{GraphFullCIs.XY}
\title{
Plot all the DR and BLS potential solutions}
\description{
Display a plot with all the DR and BLS potential solutions from OLSv to OLSh (useful for unreplicated data), choose between all the slopes (and their confidence intervals), all the intercepts (and their confidence intervals), all the confidence region (ellipses), the two extreme confidence intervals (for the expectation of Y) or the two extreme confidence bands.}
\usage{
GraphFullCIs.XY(FullCIs = NULL, BLS.estimate = NULL, lambda = NULL,
                xname = "X", yname = "Y", antilog = NULL, graph = "joint.ellipse",
                accept.int = 0, accept.int.perc = FALSE, accept.int.threshold = NULL,
                include.H0 = TRUE, include.int = TRUE)
}
\arguments{
  \item{FullCIs}{ a CIs.XY class object.}
  \item{BLS.estimate}{a BLS class object. The BLS estimates (obtained, for example, with replicated data) will be superimposed on the plot.}
  \item{lambda}{a lambdas class object. The value of \ifelse{html}{\out{&lambda;_{XY}}}{\eqn{\lambda_{XY}}{lambda_{XY}}} will be superimposed on the plot with its confidence interval.}
  \item{xname}{ a character string for the name of the X device.}
  \item{yname}{ a character string for the name of the Y device.}
  \item{antilog}{ a character string or a numeric value. This argument displays the BLS results on the initial scales on the (X,Y) plot if a logarithmic transformation was used prior to the BLS function. Options available are: 10 or "e".}
  \item{graph}{ a character string for the type of graph: "slope" to plot all the slopes, "intercept" to plot all the intercepts,"joint.ellipse" to plot all the joint confidence intervals (ellipses), "CI" to plot the confidence intervals,"CB" to plot the confidence bands, "all" to plot all the plots on different windows.}
  \item{accept.int}{ a numeric vector (length equal 1 or 2) for the value of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}}: Y = X \ifelse{html}{\out{&#177;}}{\eqn{\pm}{plus or minus}} \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} to assess whether the two devices (X and Y) are equivalent or not. Two values of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} can be entered to take into account the case where the equivalence threshold changes along the X axis.}
  \item{accept.int.perc}{ a logical variable (TRUE or FALSE) whether \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} needs to be interpreted in percentage: Y = X \ifelse{html}{\out{&#177;}}{\eqn{\pm}{plus or minus}} \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}}\%.}
  \item{accept.int.threshold}{ a numeric value for the X threshold where the value of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} changes if two \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} values are used in \code{accept.int}.}
  \item{include.H0}{a logical variable (TRUE or FALSE) whether the null hypothesis (slope = 1, intercept = 0) should lie on plot (in the case that the null hypothesis is out from the axes limits).}
  \item{include.int}{a logical variable (TRUE or FALSE) whether the confidence intervals should lie entirely in the plot.}
}
\details{The ellipses are plotted in an (\ifelse{html}{\out{&beta;}}{\eqn{\beta}{beta}},\ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}}) coordinate system where the acceptance interval is a diamond. The slopes and the intercepts are plotted on the Y-axis with \ifelse{html}{\out{&lambda;_{XY}}}{\eqn{\lambda_{XY}}{lambda_{XY}}} assigned on the X-axis. The confidence intervals and confidence bands are displayed on a classical (X,Y) plot.}
\value{ The plot requested by the argument \code{graph}.}
\references{
Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems, 2014; 134:123-139.
\cr Francq BG, Govaerts BB. Hyperbolic confidence bands of errors-in-variables regression lines applied to method comparison studies. Journal de la Societe Francaise de Statistique 2014; 155(1):23-45.}
\author{
Bernard G FRANCQ}
\seealso{
\code{\link{XY.plot}}, \code{\link{BLS}}}
\examples{
library(BivRegBLS)
data(SBP)
# Estimate all the solutions with the DR and BLS regressions
res.full=FullCIs.XY(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
# Estimate the BLS regression with replicated data
res.BLS=BLS(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
# Estimate the measurement error variances ratio
res.lambda=lambdas(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
# Plot all the results with all the potential solutions
#   and superimpose the BLS regression estimated with replicated data
GraphFullCIs.XY(FullCIs=res.full,BLS.estimate=res.BLS,lambda=res.lambda,
                xname="J",yname="S",graph="all",accept.int=10,accept.int.perc=FALSE,include.H0=TRUE)
# Plot all the potential solutions for the confidence intervals and add two acceptance intervals
GraphFullCIs.XY(FullCIs=res.full,xname="J",yname="S",graph="CI",
                accept.int=c(8,12),accept.int.threshold=150,accept.int.perc=FALSE,include.H0=TRUE)
data(Aromatics)
# Estimate all the solutions with the BLS regression
res.full=FullCIs.XY(data=Aromatics,xcol=3,ycol=4)
# Plot all the potential solutions for the confidence regions (ellipses)
GraphFullCIs.XY(FullCIs=res.full,xname="HPLC",yname="GC MS",graph="joint.ellipse")
}
\keyword{Deming Regression}
\keyword{DR}
\keyword{Bivariate Least Square regression}
\keyword{BLS}
\keyword{ellipse}
\keyword{confidence region}
