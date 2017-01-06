\name{MD.horiz.lines}
\alias{MD.horiz.lines}
\title{
Tolerance intervals in a (M,D) plot}
\description{
Calculate agreement interval and tolerance intervals (bounded by 2 horizontal lines) for the differences (D=Y-X) in a (M,D) plot (also called Bland-Altman plot).}
\usage{
MD.horiz.lines(data = NULL, xcol = 1, ycol = 2,
               pred.level = 0.95, TI.conf.level = 0.8)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X column(s) or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y column(s) or a character vector with the column names.}
  \item{pred.level}{ a numeric value for the predictive level (expressed between 0 and 1) of the tolerance intervals (beta expectation (type I) and beta gamma content (type II)).}
  \item{TI.conf.level}{ a numeric value for the confidence level of the beta gamma content tolerance interval (tolerance interval type II).}
}
\details{
The data argument is mandatory while other arguments are optional. If the data are replicated, the tolerance intervals predict where a future single difference (Di=Yi-Xi) will lie (and not an average difference). Tolerance intervals are better (than agreement interval) and should be preferred. The tolerance intervals are calculated on the univariate distribution of the differences (Di). These intervals are valid under the assumption that there is no proportional bias. If a pattern is observed, the \code{\link{CBLS}} function (CBLS regression) must be used (with its predictive intervals).}
\value{
A MD.horiz.lines class object: a list including the following elements:
  \item{data.MD}{ a table with the means ((X+Y)/2) and differences (Y-X).}
  \item{Table.Differences}{ a table with one row and several descriptive statistics: the mean of the differences, the standard deviation of the mean difference (Di=Yi-Xi) and the standard deviation for a single difference (Yik-Xik), the minimum, 1st quartile, median, 3rd quartile, maximum and number of observations.}
  \item{Intervals.horiz.lines}{ a table with 3 rows for the agreement interval, the beta expectation tolerance interval and the beta gamma content tolerance interval, with their interpretation (columns).}
}
\references{
Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.}
\author{Bernard G FRANCQ}
\seealso{
\code{\link{CBLS}}
}
\examples{
library(BivRegBLS)
data(SBP)
res.MD.horiz=MD.horiz.lines(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10,pred.level=0.95)
res.MD.horiz$Intervals.horiz.lines
}
\keyword{tolerance interval}
\keyword{MD plot}
\keyword{Bland-Altman}
\keyword{agreement}
