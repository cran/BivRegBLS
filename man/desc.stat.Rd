\name{desc.stat}
\alias{desc.stat}
\title{
Descriptive statistics in method comparison studies}
\description{
Calculate several descriptive statistics in method comparison studies per device (X and Y) and per type of samples.}
\usage{
desc.stat(data = NULL, xcol = 1, ycol = 2, IDcol = NULL)
}
\arguments{
  \item{data}{ a data set (data frame or matrix).}
  \item{xcol}{ a numeric vector to specify the X column(s) or a character vector with the column names.}
  \item{ycol}{ a numeric vector to specify the Y column(s) or a character vector with the column names.}
  \item{IDcol}{ a numeric or character variable to specify the column with the different IDs or type of samples.}
}
\details{
If \code{IDcol} is null (as by default), the descriptive statistics are calculated for X and Y. Otherwise, the descriptive statistics are calculated for X and Y for each type of sample (each ID) (with a maximum of 30 different IDs). This information is also used to differentiate the observations on a raw plot when the function \code{\link{raw.plot}} is used. In presence of missing values on X or Y and non-replicates, the rows with missing values are removed. In presence of replicates, the rows with missing values are removed if all Xi or all Yi are missing. 
\cr The results (\code{Xij}, \code{Yik}, \code{Xi}, \code{Yi}, \code{nxi}, \code{nyi}, \code{variances_x}, \code{variances_y}) are reordered according to the increasing values of Xi (the X mean values).}
\value{
A list including the following elements:
  \item{Xij}{ a table with the (replicated) X measurements (replicates are in columns).}
  \item{Yik}{ a table with the (replicated) Y measurements (replicates are in columns).}
  \item{Xi}{ a vector with the means of the X measurements.}
  \item{Yi}{ a vector with the means of the Y measurements.}
  \item{IDs}{ a vector with the different IDs.}
  \item{nxi}{ a vector with the number of X replicates per sample (patient).}
  \item{nyi}{ a vector with the number of Y replicates per sample (patient).}
  \item{variances_x}{ a vector with the variances calculated on the X replicates per sample (patient).}
  \item{variances_y}{ a vector with the variances calculated on the Y replicates per sample (patient).}
  \item{Order.Xi}{ a vector with the order of the means of the X replicates.}
  \item{statistics}{ a table with different descriptive statistics per type of sample (rows): the number of sample (patient), the number of replicates in X and Y, the degrees of freedom of the measurement error variances in X and Y, the mean, the sum of squares (Sxx and Syy), the cross-product (Sxy), the variance, minimum, 1st quartile, median, 3rd quartile, maximum for X and Y, and the Pearson correlation coefficient and its square.}
}
\references{
Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.
\cr Francq BG. Errors-in-variables regressions to assess equivalence in method comparison studies. Ph.D. Thesis, Universite Catholique de Louvain, Institute of Statistics, Biostatistics and Actuarial science, Louvain-la-Neuve, Belgium, 2013.
}
\author{Bernard G FRANCQ}
\seealso{\code{\link{lambdas}}, \code{\link{raw.plot}}
}
\examples{
library(BivRegBLS)
data(Aromatics)
res=desc.stat(data=Aromatics,xcol=3,ycol=4,IDcol=2)
}
\keyword{ descriptive statistics}
\keyword{ method comparison studies}
