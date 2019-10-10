\name{BivRegBLS-package}
\alias{BivRegBLS-package}
\alias{BivRegBLS}
\docType{package}
\title{
Tolerance Intervals and Errors-in-Variables Regressions in Method Comparison Studies (Bridging Studies)}
\description{
Assess the agreement in method comparison studies by tolerance intervals and (correlated)-errors-in-variables regressions. The results can then be displayed in an (X,Y) plot or (M,D) plot (Bland-Atman plot). The vertical and horizontal Ordinary Least Square regressions (OLSv, OLSh), the Deming Regression (DR), and the (Correlated)-Bivariate Least Square regressions (BLS and CBLS) can be used with unreplicated or replicated data. The measurement error variances and their ratio (lambda) can be estimated by an unbiased estimator with replicated data. If lambda is unknown and not estimatable, there is no unique solution and all the potential solutions can be calculated from OLSv to OLSh in a (X,Y) plot (extreme solutions), or equivalently from a correlation (between the measurement errors in a Bland-Altman plot) -1 to +1 with the CBLS in a (M,D) plot. The \code{\link{BLS}} and \code{\link{CBLS}} are the two main regressions. They provide a table of the estimates (estimates, standard error, confidence intervals and pvalues for separate and joint hypotheses), the coordinates of the joint confidence interval (confidence region, or ellipse), and the four following hyperbolic intervals: the Confidence Intervals (CI), the Confidence Bands (CB), the Predictive Intervals (PI), and the Generalized predictive Intervals (GI). The \code{\link{XY.plot}} and \code{\link{MD.plot}} are the two main graphical functions to display an (X,Y) plot or (M,D) plot with the desired results. If one can assume no proportional bias, the (M,D) plot may be simplified by calculating horizontal lines intervals with the beta-expectation tolerance interval (type I) or the beta-gamma content tolerance interval (type II).
}
\details{
\tabular{ll}{
Package: \tab BivRegBLS\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2019-10-10\cr
License: \tab AGPL-3\cr
}
The most important functions are \code{\link{BLS}} (Bivariate Least Square regression) and \code{\link{CBLS}} (Correlated Bivariate Least Square regression). The results can then be plotted with respectively the functions \code{\link{XY.plot}} and \code{\link{MD.plot}}. Univariate tolerance intervals (bounded by two straight lines in the (M,D) plot) for the differences by two measurement methods can be obtained by the function \code{\link{MD.horiz.lines}}.}
\author{
Bernard G FRANCQ <BivRegBLS@gmail.com>, Marion BERGER <marion.berger@sanofi.com>
\cr Maintainer: \packageMaintainer{BivRegBLS}
}
\note{
BivRegBLS was developed with a partnership between the University of Glasgow and Sanofi under a Knowledge Exchange award: BGF lead researcher, MB project manager, and big thanks to the whole project team during this fantastic and wonderful journey: Christophe Agut, Armand Berges, Guy Mathieu, Franck Pellissier, Veronique Onado and Delphine Attonaty.
}
\references{
Francq BG, Govaerts BB. How to regress and predict in a Bland-Altman plot? Review and contribution based on tolerance intervals and correlated-errors-in-variables models. Statistics in Medicine, 2016; 35:2328-2358.
\cr Francq BG, Govaerts BB. Measurement methods comparison with errors-in-variables regressions. From horizontal to vertical OLS regression, review and new perspectives. Chemometrics and Intelligent Laboratory Systems, 2014; 134:123-139.
\cr Francq BG, Govaerts BB. Hyperbolic confidence bands of errors-in-variables regression lines applied to method comparison studies. Journal de la Societe Francaise de Statistique 2014; 155(1):23-45.
\cr Francq BG. Errors-in-variables regressions to assess equivalence in method comparison studies. Ph.D. Thesis, UCLouvain, Institute of Statistics, Biostatistics and Actuarial science, Louvain-la-Neuve, Belgium, 2013.
}
\keyword{BivRegBLS}
\examples{
library(BivRegBLS)
data(SBP)
### Descriptive statistics
res=desc.stat(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
raw.plot(data.plot=res,xname="J",yname="S",graph="XY.bar.SEM")
### BLS regression in an (X,Y) plot
res.BLS=BLS(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10,qx=3,qy=3)
XY.plot(BLS.results=res.BLS,xname="J",yname="S",accept.int=10,accept.int.perc=FALSE)
### CBLS regression in an (M,D) plot
res.CBLS=CBLS(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10)
MD.plot(results=res.CBLS,xname="J",yname="S",
        accept.int=10,accept.int.perc=FALSE,include.int=TRUE,graph.int=c("CI","GI"))
### Univariate tolerance intervals in an (M,D) plot
res.MD.horiz=MD.horiz.lines(data=SBP,xcol=c("J1","J2","J3"),ycol=8:10,pred.level=0.95)
MD.plot(results=res.MD.horiz,xname="J",yname="S",
        accept.int=10,accept.int.perc=FALSE,include.int=TRUE,graph.horiz.int=c("bTI","bgTI"))
}
