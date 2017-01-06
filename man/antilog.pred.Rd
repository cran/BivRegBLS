\name{antilog.pred}
\alias{antilog.pred}
\title{Back transforms the results if a logarithmic transformation is used}
\description{
If the data are log-normal, the user can apply a logarithmic transformation. Then, \code{\link{antilog.pred}} will automatically back-transform (exponentiate) the data and the predictions (hyperbolic intervals) obtained by the \code{\link{BLS}}, \code{\link{CBLS}}, \code{\link{MD.horiz.lines}}, \code{\link{FullCIs.XY}} or \code{\link{FullCIs.MD}} functions.}
\usage{
antilog.pred(results = NULL, base = 10)
}
\arguments{
  \item{results}{an object of class \code{BLS}, \code{CBLS}, \code{MD.horiz.lines}, \code{CIs.XY} or \code{CIs.MD}.}
  \item{base}{ a character string or a numeric value. Options available are: 10 (default value) or "e".}
}
\details{\code{\link{antilog.pred}} is automatically called by the argument \code{antilog} in the functions \code{\link{XY.plot}} or \code{\link{MD.plot}}.}
\value{An object of class \code{BLS}, \code{CBLS}, \code{MD.horiz.lines}, \code{CIs.XY} or \code{CIs.MD}.}

\author{Bernard G FRANCQ}
\seealso{
\code{\link{BLS}}, \code{\link{CBLS}}, \code{\link{MD.horiz.lines}}, \code{\link{FullCIs.XY}}, \code{\link{FullCIs.MD}}}
\examples{
library(BivRegBLS)
data(SBP)
SBPlog=SBP
SBPlog[,2:10]=log(SBPlog[,2:10])
res.BLS.log=BLS(data=SBPlog,xcol=c("J1","J2","J3"),ycol=8:10)
res.BLS=antilog.pred(results=res.BLS.log,base="e")
}
\keyword{log normal}
\keyword{exponential}
