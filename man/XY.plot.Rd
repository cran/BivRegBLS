\name{XY.plot}
\alias{XY.plot}
\title{
Display the BLS regression in a (X,Y) plot
}
\description{
Display the BLS regression in a (X,Y) plot with or without hyperbolic confidence and/or predictive intervals, and an acceptance interval.}
\usage{
XY.plot(BLS.results = NULL, xname = "X", yname = "Y", antilog = NULL,
        accept.int = 0, accept.int.perc = FALSE, accept.int.threshold = NULL,
        graph.int = c("CB", "PI"), include.int = FALSE, col.BLS = 1, col.CI = 2,
        col.CB = 3, col.PI = 4, col.GI = 5, lty.BLS = 1, lty.CI = 1,
        lty.CB = 1, lty.PI = 1, lty.GI = 1, ...)
}
\arguments{
  \item{BLS.results}{a BLS or BLS.ht class object (obtained with the \code{\link{BLS}} or \code{\link{BLS.ht}} functions).}
  \item{xname}{ a character string for the name of the X device.}
  \item{yname}{ a character string for the name of the Y device.}
  \item{antilog}{ a character string or a numeric value. This argument display the BLS results on the initial scales on the (X,Y) plot if a logarithmic transformation was used prior to the BLS function. Options available are: 10 or "e".}
  \item{accept.int}{ a numeric vector (length equal 1 or 2) for the value of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}}: Y = X \ifelse{html}{\out{&#177;}}{\eqn{\pm}{plus or minus}} \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} to assess whether two devices (X and Y) are equivalent or not. Two values of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} can be entered to take into account the case where the equivalence threshold changes along the X axis.}
  \item{accept.int.perc}{ a logical variable (TRUE or FALSE) whether \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} needs to be interpreted in percentage: Y = X \ifelse{html}{\out{&#177;}}{\eqn{\pm}{plus or minus}} \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} \%.}
  \item{accept.int.threshold}{ a numeric value for the X threshold where the value of \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} changes if two \ifelse{html}{\out{&Delta;}}{\eqn{\Delta}{Delta}} values are used in \code{accept.int}.}
  \item{graph.int}{ a character vector for the hyperbolic curves to be displayed on the graph. Options are: "CI" (Confidence Intervals),"CB" (Confidence Bands),"PI" (Predictive Intervals),"GI" (Generalised Intervals).}
  \item{include.int}{ a logical variable whether the hyperbolic intervals should lie entirely in the plot.}
  \item{col.BLS}{ a character string or a numeric variable for the colour of the BLS regression line.}
  \item{col.CI}{ a character string or a numeric variable for the colour of the confidence intervals.}
  \item{col.CB}{ a character string or a numeric variable for the colour of the confidence bands.}
  \item{col.PI}{ a character string or a numeric variable for the colour of the predictive intervals.}
  \item{col.GI}{ a character string or a numeric variable for the colour of the generalized intervals.}
  \item{lty.BLS}{ a numeric variable for the type of line of the BLS regression line.}
  \item{lty.CI}{ a numeric variable for the type of line for the confidence intervals.}
  \item{lty.CB}{ a numeric variable for the type of line for the confidence bands.}
  \item{lty.PI}{ a numeric variable for the type of line for the predictive intervals.}
  \item{lty.GI}{ a numeric variable for the type of line for the generalized intervals.}
  \item{\dots}{ the common arguments from 'plot' or 'par' that may be used additionaly, such as xlim, ylim, xlab, ylab.}
}
\details{The BLS.result argument is mandatory.}
\value{An (X,Y) plot in a new window.}
\references{
Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.
\cr Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems, 2014; 134:123-139.
\cr Francq BG, Govaerts BB. Hyperbolic confidence bands of errors-in-variables regression lines applied to method comparison studies. Journal de la Societe Francaise de Statistique 2014; 155(1):23-45.}
\author{Bernard G FRANCQ}
\note{
The limits of the axes and their labels are set automatically by the function. To compare different plots with fixed limits, use xlim and ylim. To write customized labels, use xlab and ylab.}
\seealso{
\code{\link{MD.plot}}
}
\examples{
library(BivRegBLS)
data(SBP)
# Estimate the BLS regression on replicated data
res.BLS=BLS(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10,qx=3,qy=3)
# Plot the results in a (X,Y) plot with an acceptance interval
XY.plot(BLS.results=res.BLS,xname="J",yname="S",accept.int=10,accept.int.perc=FALSE)
XY.plot(BLS.results=res.BLS,xname="J",yname="S",accept.int=10,accept.int.perc=TRUE)
XY.plot(BLS.results=res.BLS,xname="J",yname="S",accept.int=c(10,20),
        accept.int.perc=FALSE,accept.int.threshold=150)
}
\keyword{BLS regression}
\keyword{XY plot}
\keyword{confidence interval}
\keyword{predictive interval}
\keyword{confidence bands}
\keyword{generalized interval}
