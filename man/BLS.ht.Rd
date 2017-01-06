\name{BLS.ht}
\alias{BLS.ht}
\title{
Bivariate Least Square regression (BLS)}
\description{
Estimate the heteroscedastic Bivariate Least Square regression with replicated data or variance functions.}
\usage{
BLS.ht(data = NULL, xcol = NULL, ycol = NULL, var.x.col = NULL, var.y.col = NULL,
       var.x.formula = NULL, var.y.formula = NULL, nxi.col = NULL, nyi.col = NULL,
       df.var.x.col = NULL, df.var.y.col = NULL, conf.level = 0.95, pred.level = 0.95,
       npoints = 1000, qx = 1, qy = 1, xpred = NULL, tol = 1e-05)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X column(s) or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y column(s) or a character vector with the column names.}
  \item{var.x.col}{ a numeric or character variable to specify the column with the heteroscedastic variances of X.}
  \item{var.y.col}{ a numeric or character variable to specify the column with the heteroscedastic variances of Y.}
  \item{var.x.formula}{ a character string with a formula related to X that will be interpreted as a variance function, i.e. var.x.formula="exp(2*X+3)".}
  \item{var.y.formula}{ a character string with a formula related to Y that will be interpreted as a variance function, i.e. var.y.formula="exp(2*Y+3)".}
  \item{nxi.col}{ a numeric or character variable to specify the column with the number of replicates in X on which the variances were calculated.}
  \item{nyi.col}{ a numeric or character variable to specify the column with the number of replicates in Y on which the variances were calculated.}
  \item{df.var.x.col}{ a numeric or character variable to specify the column with the degrees of freedom of the heteroscedastic variances in X.}
  \item{df.var.y.col}{ a numeric or character variable to specify the column with the degrees of freedom of the heteroscedastic variances in Y.}
  \item{conf.level}{ a numeric value for the confidence level (expressed between 0 and 1).}
  \item{pred.level}{ a numeric value for the predictive level (expressed between 0 and 1).}
  \item{npoints}{ a numeric variable for the number of points to smooth the hyperbolic curves.}
  \item{qx}{ a numeric value to predict the mean of qy Y future values from the mean of qxX values (generalized interval).}
  \item{qy}{ a numeric value to predict the mean of qy Y future values from the mean of qx X values (generalized interval).}
  \item{xpred}{ a numeric vector for customized predictions at given X values.}
  \item{tol}{ a numeric variable to change the tolerance for the BLS to converge.}
}
\details{
The data argument is mandatory. If the data are unreplicated, then the measurement error variances must be given. The confidence level is used for the confidence intervals of the parameters (\ifelse{html}{\out{&beta;}}{\eqn{\beta}{beta}}, \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}}), the hyperbolic confidence intervals  (the prediction of the expectation of Y for a given X) and the hyperbolic confidence bands. The predictive level is used for the hyperbolic predictive intervals (the prediction of a future Y for a given X) and the hyperbolic generalized intervals (the prediction of the mean of qy future Y values from a given (mean of) X).
\cr The results (\code{Xij}, \code{Yik}, \code{Xi}, \code{Yi}, \code{nxi}, \code{nyi}, \code{variances_x}, \code{variances_y}) are reordered according to the increasing values of Xi (the X mean values).
}
\value{
A BLS.ht class object, a list including the following elements:

  \item{Xij}{ a table with the (replicated) X measurements (replicates are in columns).}
  \item{Yik}{ a table with the (replicated) Y measurements (replicates are in columns).}
  \item{Xi}{ a vector with the means of the X measurements.}
  \item{Yi}{ a vector with the means of the Y measurements.}
  \item{nxi}{ a vector with the number of X replicates per sample (patient).}
  \item{nyi}{ a vector with the number of Y replicates per sample (patient).}
  \item{variances_x}{ a vector with the variances calculated on the X replicates per sample (patient).}
  \item{variances_y}{ a vector with the variances calculated on the Y replicates per sample (patient).}
  \item{Ellipse.BLS.ht}{ a two columns matrix with the coordinates of the joint confidence interval  (confidence region, ellipse) for the parameters (\ifelse{html}{\out{&beta;}}{\eqn{\beta}{beta}}, \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}}).}
  \item{Estimate.BLS.ht}{ a table (data frame) with the estimates of the intercept and the slope, standard error, confidence interval and pvalue (null hypothesis: slope = 1, intercept = 0).}
  \item{Pred.BLS.ht}{ a data frame with npoints rows (from the minimum to the maximum of the observed X values) and the following columns: the X values where the predictions are calculated (X0), the Y predicted values (Ypred), the lower and upper bounds of the confidence interval, predictive interval, generalized interval and confidence bands.}
  \item{xpred.BLS.ht}{ a data frame with the customized predictions and the same columns than \code{Pred.BLS.ht}.}
}
\note{
The prediction interval should be interpreted with caution as it is still under development.
}
\references{
Francq BG. Errors-in-variables regressions to assess equivalence in method comparison studies. Ph.D. Thesis, Universite Catholique de Louvain, Institute of Statistics, Biostatistics and Actuarial science, Louvain-la-Neuve, Belgium, 2013.}
\author{Bernard G FRANCQ}

\seealso{
\code{\link{BLS}}
}
\examples{
library(BivRegBLS)
data(SBP)
res.BLS.ht=BLS.ht(data=SBP,xcol=c("J1","J2","J3"),ycol=c("S1","S2","S3"))
}
\keyword{BLS}
\keyword{Bivariate Least Square}
\keyword{heteroscedastic}
\keyword{regression}
\keyword{confidence region}
\keyword{joint confidence interval}
\keyword{ellipse}

