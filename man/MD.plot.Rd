\name{MD.plot}
\alias{MD.plot}
\title{
Display the CBLS regression, or univariate tolerance intervals in a (M,D) plot}
\description{
Display the CBLS regression in a (M,D) plot with or without hyperbolic confidence and/or predictive intervals, and an acceptance interval. Alternatively, univariate tolerance intervals which are bounded by two horizontal lines can be plotted.}
\usage{
MD.plot(results = NULL, xname = "X", yname = "Y", antilog = NULL,
        accept.int = 0, accept.int.perc = FALSE, accept.int.threshold = NULL,
        include.int = TRUE, graph.int = c("CB", "PI"), graph.horiz.int = c("bTI", "bgTI"),
        col.CBLS = 1, col.CI = 2, col.CB = 3, col.PI = 4, col.GI = 5, col.bTI = 3,
        col.bgTI = 4, lty.CBLS = 1, lty.CI = 1, lty.CB = 1, lty.PI = 1, lty.GI = 1,
        lty.bTI = 1, lty.bgTI = 1, ...)
}
\arguments{
  \item{results}{a CBLS class object (obtained with the \code{\link{BLS}} function) or a MD.horiz.lines class object (obtained with the \code{\link{MD.horiz.lines}} function).}
  \item{xname}{ a character string for the name of the X device.}
  \item{yname}{ a character string for the name of the Y device.}
  \item{antilog}{ a character string or a numeric value. This argument displays the CBLS or MD.horiz.lines results on the initial scales in the (M,D) plot if a logarithmic transformation was used prior to the CBLS or MD.horiz.lines functions. Options available are: 10 or "e".}
  \item{accept.int}{ a numeric vector (length equal 1 or 2) for the value of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}}: |Y-X|<\ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} to assess whether two devices (X and Y) are equivalent or not. Two values of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} can be entered to take into account the case where the equivalence threshold changes along the M axis.}
  \item{accept.int.perc}{ a logical variable (True or Flase) whether \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} needs to be interpreted in percentage: Y = X \ifelse{html}{\out{&#177;}}{\eqn{\pm}{plus or minus}} \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}}\%.}
  \item{accept.int.threshold}{ a numeric value for the M threshold where the value of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} changes if two \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} values are used in \code{accept.int}.}
  \item{include.int}{ a logical variable (True or Flase) whether the hyperbolic intervals should lie entirely in the plot.}
  \item{graph.int}{ a character vector for the hyperbolic curves to be displayed on the graph. Options are: "CI" (Confidence Intervals),"CB" (Confidence Bands),"PI" (Predictive Intervals),"GI" (Generalised Intervals).}
  \item{graph.horiz.int}{ if a MD.horiz.lines object is used, a character vector for the univariate tolerance intervals to be displayed on the graph. Options are: "bTI" (beta expectation tolerance interval or tolerance interval type I) and/or "bgTI" (beta gamma content tolerance interval or tolerance interval type II).}
  \item{col.CBLS}{ a character string or a numeric variable for the colour of the CBLS regression line.}
  \item{col.CI}{ a character string or a numeric variable for the colour of the confidence intervals.}
  \item{col.CB}{ a character string or a numeric variable for the colour of the confidence bands.}
  \item{col.PI}{ a character string or a numeric variable for the colour of the predictive intervals.}
  \item{col.GI}{ a character string or a numeric variable for the colour of the generalized intervals.}
  \item{col.bTI}{ a character string or a numeric variable for the colour of the beta expectation tolerance interval.}
  \item{col.bgTI}{ a character string or a numeric variable for the colour of the beta gamma content tolerance interval.}
  \item{lty.CBLS}{ a numeric variable for the type of line of the CBLS regression line.}
  \item{lty.CI}{ a numeric variable for the type of line for the confidence intervals.}
  \item{lty.CB}{ a numeric variable for the type of line for the confidence bands.}
  \item{lty.PI}{ a numeric variable for the type of line for the predictive intervals.}
  \item{lty.GI}{ a numeric variable for the type of line for the generalized intervals.}
  \item{lty.bTI}{ a numeric variable for the type of line for the beta expectation tolerance interval.}
  \item{lty.bgTI}{ a numeric variable for the type of line for the beta gamma content tolerance interval.}
  \item{\dots}{ the common arguments from 'plot' or 'par' that may be used additionaly, such as xlim, ylim, xlab, ylab.}
}
\details{
The \code{results} argument is mandatory. The value of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} (\code{accept.int}) is converted to percentage if antilog is used in a (M,D) plot to define 2 asymetric bounds (1- \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}}/100, (100/(100-\ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}})).}
\value{An (M,D) plot in a new window.}
\references{
Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.
\cr Francq BG, Govaerts BB. Hyperbolic confidence bands of errors-in-variables regression lines applied to method comparison studies. Journal de la Societe Francaise de Statistique 2014; 155(1):23-45.}
\author{Bernard G FRANCQ}
\note{The limits of the axes and their labels are set automatically by the function. To compare different plots with fixed limits, use xlim and ylim. To write customized labels, use xlab and ylab.}
\seealso{
\code{\link{XY.plot}}
}
\examples{
library(BivRegBLS)
data(SBP)
# Estimate the CBLS regression on replicated data
res.CBLS=CBLS(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
# Plot the results in a (M,D) plot with an acceptance interval
MD.plot(results=res.CBLS,xname="J",yname="S",accept.int=10,accept.int.perc=FALSE)
MD.plot(results=res.CBLS,xname="J",yname="S",accept.int=10,accept.int.perc=TRUE)
MD.plot(results=res.CBLS,xname="J",yname="S",accept.int=c(10,15),
        accept.int.perc=FALSE,accept.int.threshold=150)
}
\keyword{CBLS regression}
\keyword{MD plot}
\keyword{Bland-Altman}
\keyword{confidence interval}
\keyword{predictive interval}
\keyword{confidence bands}
\keyword{generalized interval}
