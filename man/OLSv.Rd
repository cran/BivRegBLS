\name{OLSv}
\alias{OLSv}
\title{Vertical Ordinary Least Square regression}
\description{
Fit a linear ordinary least square regression by minimising the residuals in a vertical direction.}
\usage{
OLSv(data = NULL, xcol = 1, ycol = 2, conf.level = 0.95, pred.level = 0.95,
     npoints = 1000, q = 1, xpred = NULL)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X column(s) or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y column(s) or a character vector with the column names.}
  \item{conf.level}{ a numeric value for the confidence level (expressed between 0 and 1).}
  \item{pred.level}{ a numeric value for the predictive level (expressed between 0 and 1).}
  \item{npoints}{ an integer (at least 10) for the number of points to smooth the hyperbolic curves.}
  \item{q}{ an integer to predict the mean of q future values (generalized interval).}
  \item{xpred}{ a numeric vector for customized predictions at given X values.}
}
\details{
The data argument is mandatory while other arguments are optional. The confidence level is used for the confidence intervals of the parameters, the hyperbolic confidence intervals  (the prediction of the expectation of Y for a given X) and the hyperbolic confidence bands. The predictive level is used for the hyperbolic predictive intervals (the prediction of a future Y for a given X) and the hyperbolic generalized intervals (the prediction of the mean of q future Y values for a given X).
}
\value{
A list including the following elements:
  \item{Ellipse.OLSv}{ a two columns matrix with the coordinates of the joint confidence interval (confidence region, ellipse) for the parameters (\ifelse{html}{\out{&beta;}}{\eqn{\beta}{beta}}, \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}}).}
  \item{Estimate.OLSv}{ a table (data frame) with the estimates of the intercept and the slope, standard error, confidence interval and pvalue (null hypothesis: slope = 1, intercept = 0).}
  \item{Pred.OLSv}{ a data frame with npoints rows (from the minimum to the maximum of the observed X values) and the following columns: the X values where the predictions are calculated (X0), the Y predicted values (Ypred), the lower and upper bounds of the confidence intervals, predictive intervals, generalized intervals and confidence bands.}
  \item{xpred.OLSv}{ a data frame with the customized predictions and the same columns than \code{Pred.OLSv}.}
}
\references{
Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems 2014; 134:123-139.
\cr Francq BG, Govaerts BB. Hyperbolic confidence bands of errors-in-variables regression lines applied to method comparison studies. Journal de la Societe Francaise de Statistique 2014; 155(1):23-45.
}
\author{Bernard G FRANCQ}
\note{ The default value for xcol (ycol) is 1 (2) for the 1st (2nd) column.}
\seealso{
\code{\link{OLSh}}
}
\examples{
res.OLSv=OLSv(matrix(nrow=10,ncol=2,c(1:10,(1:10)+rnorm(10))))
res.OLSv$Estimate.OLSv
}

\keyword{regression}
\keyword{ordinary least square}
\keyword{vertical}
\keyword{confidence interval}
\keyword{confidence region}
\keyword{ellipse}
\keyword{confidence bands}
\keyword{predictive interval}