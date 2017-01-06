\name{CBLS}
\alias{CBLS}
\title{
Correlated Bivariate Least Square regression (CBLS)}
\description{
Estimate the Correlated Bivariate Least Square regression with replicated data (in a (M,D) plot) where M=(X+Y)/2 and D=Y-X.}
\usage{
CBLS(data = NULL, xcol = 1, ycol = 2, var.x = NULL, var.y = NULL,
     df.var.x = Inf, df.var.y = Inf, ratio.var = NULL, conf.level = 0.95,
     pred.level = 0.95, npoints = 1000, qx = 1, qy = 1, xpred = NULL)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X column(s) or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y column(s) or a character vector with the column names.}
  \item{var.x}{ a numeric variable for the variance of the measurement error of device X if known.}
  \item{var.y}{ a numeric variable for the variance of the measurement error of device Y if known.}
  \item{df.var.x}{ a numeric variable for the degrees of freedom of the variance of the measurement error of device X if known.}
  \item{df.var.y}{ a numeric variable for the degrees of freedom of the variance of the measurement error of device Y if known.}
  \item{ratio.var}{ a numeric value for \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}, the ratio of the measurement error variances (Y over X) if known.}
  \item{conf.level}{ a numeric value for the confidence level (expressed between 0 and 1).}
  \item{pred.level}{ a numeric value for the predictive level (expressed between 0 and 1).}
  \item{npoints}{ a numeric variable for the number of points to smooth the hyperbolic curves.}
  \item{qx}{ a numeric value to predict the mean of qy Y future values from the mean of qx X values (generalized interval).}
  \item{qy}{ a numeric value to predict the mean of qy Y future values from the mean of qx X values (generalized interval).}
  \item{xpred}{ a numeric vector for customized predictions at given M values.}
}
\details{
The data argument is mandatory. If the data are unreplicated, then the measurement error variances must be given or their ratio (\ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}) in order to calculate the correlation, \ifelse{html}{\out{&rho;_{MD}}}{\eqn{\rho_{MD}}{rho_{MD}}}, between the measurement errors of the differences (on the Y-axis) and the measurement errors of the means (on the X-axis). The confidence level is used for the confidence intervals of the parameters (\ifelse{html}{\out{&rho;_{MD}}}{\eqn{\rho_{MD}}{rho_{MD}}}, \ifelse{html}{\out{&beta;}}{\eqn{\beta}{beta}} (slope), \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} (intercept)), the hyperbolic confidence intervals  (the prediction of the expectation of Y for a given X) and the hyperbolic confidence bands. The predictive level is used for the hyperbolic predictive intervals (the prediction of a future Y for a given X) and the hyperbolic generalized intervals (the prediction of the mean of q future Y values for a given X).
}
\value{
A CBLS class object, a list including the following elements:

  \item{Xij}{ a table with the (replicated) X measurements (replicates are in columns).}
  \item{Yik}{ a table with the (replicated) Y measurements (replicates are in columns).}
  \item{Xi}{ a vector with the means of the X measurements.}
  \item{Yi}{ a vector with the means of the Y measurements.}
  \item{Mi}{ a vector with the means ((X+Y)/2) measurements.}
  \item{Di}{ a vector with the differences (Y-X) measurements.}
  \item{nxi}{ a vector with the number of X replicates per sample (patient).}
  \item{nyi}{ a vector with the number of Y replicates per sample (patient).}
  \item{variances_x}{ a vector with the variances calculated on the Y replicates per sample (patient).}
  \item{variances_y}{ a vector with the variances calculated on the Y replicates per sample (patient).}
  \item{Rho.MD}{ a table with the value of \ifelse{html}{\out{&rho;MD}}{\eqn{\rho_{MD}}{rhoXY}} (the correlation between the measurement errors of the means and the differences) and its confidence interval.}
  \item{Ellipse.CBLS}{ a two columns matrix with the coordinates of the joint confidence interval (confidence region, ellipse) for the parameters (\ifelse{html}{\out{&beta;}}{\eqn{\beta}{beta}}, \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}}).}
  \item{Estimate.CBLS}{ a table (data frame) with the estimates of the intercept and the slope, standard error, confidence interval and pvalue (null hypothesis: slope = 0, intercept = 0).}
  \item{Pred.CBLS}{ a data frame with npoints rows (from the minimum to the maximum of the observed X values) and the following columns: the X values where the predictions are calculated (X0), the Y predicted values (Ypred), the lower and upper bounds of the confidence interval, predictive interval, generalized interval and confidence bands.}
  \item{xpred.CBLS}{ a data frame with the customized predictions and the same columns than \code{Pred.CBLS}.}
}
\references{
Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.}
\author{Bernard G FRANCQ}

\seealso{
\code{\link{BLS}}
}
\examples{
library(BivRegBLS)
data(SBP)
# CBLS regression on replicated data
res.CBLS1=CBLS(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10,qx=3,qy=3,xpred=c(100,120,140,160))
# CBLS regression on unreplicated data with measurement error variances previously estimated
res.CBLS2=CBLS(data=SBP,xcol=c("J1"),ycol="S1",var.x=80,var.y=50,df.var.x=100,df.var.y=100)
}
\keyword{CBLS}
\keyword{Correlated Bivariate Least Square}
\keyword{regression}
\keyword{confidence interval}
\keyword{confidence region}
\keyword{ellipse}
\keyword{confidence bands}
\keyword{predictive interval}
