\name{OLSh}
\alias{OLSh}
\title{Horizontal Ordinary Least Square regression}
\description{
Fit a linear ordinary least square regression by minimising the residuals in a horizontal direction.}
\usage{
OLSh(data = NULL, xcol = 1, ycol = 2, conf.level = 0.95)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X columns or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y columns or a character vector with the column names.}
  \item{conf.level}{ a numeric value for the confidence level (expressed between 0 and 1).}
}
\details{
The data argument is mandatory while other arugments are optional.
}
\value{
A list including the following elements:
  \item{Ellipse.OLSh}{ a two columns matrix with the coordinates of the joint confidence interval (confidence region) for the parameters (\ifelse{html}{\out{&beta;}}{\eqn{\beta}{beta}}, \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}}).}
  \item{Estimate.OLSh}{ a table (data frame) with the estimates of the intercept and the slope, standard error, confidence interval and pvalue (null hypothesis: slope = 1, intercept = 0).}
}
\references{
Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems 2014; 134:123-139.
\cr Francq BG, Govaerts BB. Hyperbolic confidence bands of errors-in-variables regression lines applied to method comparison studies. Journal de la Societe Francaise de Statistique 2014; 155(1):23-45.
}
\author{
Bernard G FRANCQ
}
\note{
 The default value for xcol (ycol) is 1 (2) for the 1st (2nd) column. The confidence region for the OLSh parameters is 'distorted' as it results from the OLSv confidence region (ellipse).
}

\seealso{
\code{\link{OLSv}}
}
\examples{
res.OLSh=OLSh(matrix(nrow=10,ncol=2,c((1:10)+rnorm(10),1:10)))
res.OLSh$Estimate.OLSh
}
\keyword{regression}
\keyword{ordinary least square}
\keyword{horizontal}
\keyword{confidence interval}
\keyword{confidence region}
\keyword{confidence bands}
