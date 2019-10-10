\name{raw.plot}
\alias{raw.plot}
\title{
Raw plot for descriptive statistics}
\description{
Display a plot with the raw data in an (X,Y) plot or (M,D) plot with or witout error bars.}
\usage{
raw.plot(data.plot = NULL, xname = "X", yname = "Y", graph = "XY.means",
         col.ID = NULL, pch.ID = NULL, ...)
}
\arguments{
  \item{data.plot}{ a desc.stat class object (obtained by using the desc.stat function).}
  \item{xname}{ a character string for the name of the X device.}
  \item{yname}{ a character string for the name of the Y device.}
  \item{graph}{ a character string for the kind of graph. Options available are: XY.means, XY.points, XY.bar.range, XY.bar.SEM, XY.bar.SD, MD.means.}
  \item{col.ID}{ a numeric or character vector describing the color of the points per ID.}
  \item{pch.ID}{ a numeric vector describing the type of points per ID to plot the data.}
  \item{\dots}{ the common arguments from 'plot' or 'par' that may be used additionaly, such as xlim, ylim, xlab, ylab.}
}
\details{
The data.plot argument is mandatory. The labels of the X and Y axes are built by default with \code{xname} and \code{yname}, and with the type of \code{graph}. The arguments \code{col.ID} and \code{pch.ID} are useful if the argument \code{IDcol} was used in the function \code{\link{desc.stat}}. \code{col.ID} and \code{pch.ID} are recycled if shorter than the number of IDs, i.e. if 3 colors are specified as col.ID=c(1,2,3) while there are 7 IDs, then the colors will be used as c(1,2,3,1,2,3,1).
\cr With the argument \code{graph}, the option \code{XY.means} plots the mean measures in a (X,Y) plot, \code{XY.points} plots all the X (Y) replicated data centered on the means of Y (X), \code{XY.bar.range} plots the error bars from the minimum to the maximum (the range) of the replicates, \code{XY.bar.SEM} plots the error bars with the standard errors of the means, \code{XY.bar.SD} plots the error bars with the standard deviations, \code{MD.means} plots the mean measures in a (M,D) plot.
}
\value{
A plot in a new window.
}
\references{
Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems 2014; 134:123-139.}
\author{
Bernard G FRANCQ}
\note{
The limits of the axes and the labels are set automatically by the function. To compare different plots with fixed limits, use xlim and ylim. To write customized labels, use xlab and ylab.}
\seealso{
\code{\link{XY.plot}}, \code{\link{MD.plot}}}
\examples{
library(BivRegBLS)
data(Aromatics)
# Calculate the descriptive statistics
res=desc.stat(data=Aromatics,xcol=3,ycol=4,IDcol="Type")
# Plot the mean or single measures (it is the same for unreplicated data)
raw.plot(data.plot=res,xname="HPLC",yname="GC MS",graph="XY.means")
raw.plot(data.plot=res,xname="HPLC",yname="GC MS",graph="XY.points")
# Plot with customized colours and type of points per type of samples
raw.plot(data.plot=res,xname="HPLC",yname="GC MS",graph="XY.points",pch.ID=c(19,5,8),col.ID=c(1,2))
raw.plot(data.plot=res,xname="HPLC",yname="GC MS",graph="MD.means")
raw.plot(data.plot=res,xname="HPLC",yname="GC MS",graph="MD.means",col.ID=c(1,2,4))
data(SBP)
# Calculate the descriptive statistics
res=desc.stat(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
# Different plots to show the measurement uncertainties
raw.plot(data.plot=res,xname="J",yname="S",graph="XY.means")
raw.plot(data.plot=res,xname="J",yname="S",graph="XY.points")
raw.plot(data.plot=res,xname="J",yname="S",graph="XY.bar.range")
raw.plot(data.plot=res,xname="J",yname="S",graph="XY.bar.SD")
raw.plot(data.plot=res,xname="J",yname="S",graph="XY.bar.SEM")
raw.plot(data.plot=res,xname="J",yname="S",graph="MD.means")
}
\keyword{plot}
\keyword{errors bars}
\keyword{XY plot}
\keyword{MD plot}
\keyword{Bland-Altman}
