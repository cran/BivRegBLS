\name{df.WS}
\alias{df.WS}
\title{
Degrees of freedom by Welch-Satterthwaite
}
\description{
Calculate the degrees of freedom from the Welch-Satterthwaite equation for a linear combination of sample variances.
}
\usage{
df.WS(variances = NULL, k = rep(1, length(variances)), dfs = NULL)
}
\arguments{
  \item{variances}{ a numeric vector for the variances.}
  \item{k}{ a numeric vector with the multiplicative constants.}
  \item{dfs}{ a numeric vector with the degrees of freedom of each variance.}
}
\details{
The \code{variances} argument is mandatory while other arguments are optional.
}
\value{
A numeric variable with the degrees of freedom of the linear combination of the variances.
}
\references{
Satterthwaite FE. An Approximate Distribution of Estimates of Variance Components. Biometrics Bulletin, 1946, 2: 110-114.
\cr Welch BK. The generalization of "student's" problem when several different population variances are involved. Biometrika, 1947, 34: 28-35.
}
\author{Bernard G FRANCQ}
\examples{
df.WS(variances=c(10,15,20),k=c(1.5,2,1.3),dfs=c(8,13,11))
}
\keyword{degrees of freedom}
\keyword{variances}
\keyword{linear combination}
\keyword{Welch-Satterwaithe}