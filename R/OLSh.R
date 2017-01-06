#' @export
OLSh <-
function(data=NULL,xcol=1,ycol=2,conf.level=0.95)
{

if(length(conf.level) > 1) stop("Only one confidence level is allowed.")
if(conf.level <= 0 | conf.level >= 1) stop("The confidence level must be between 0 and 1 (excluded)")

res=desc.stat(data=data,xcol=xcol,ycol=ycol)
x=res$Xi
y=res$Yi
xmean=res$statistics$Xmean
ymean=res$statistics$Ymean
Sxx=res$statistics$Sxx
Syy=res$statistics$Syy
Sxy=res$statistics$Sxy
n=res$statistics$N

slope_OLSh=Sxy/Syy
intercept_OLSh=xmean-ymean*slope_OLSh
S2_OLSh=sum((x-intercept_OLSh-slope_OLSh*y)^2)/(n-2)
S_slope_OLSh=sqrt(S2_OLSh/Syy)
S_intercept_OLSh=sqrt(S2_OLSh*(1/n+ymean^2/Syy))
cov_slope_intercept_OLSh=-S2_OLSh*ymean/Syy
cov_matrix_OLSh=matrix(nrow=2,ncol=2,c(S_slope_OLSh^2,cov_slope_intercept_OLSh,cov_slope_intercept_OLSh,S_intercept_OLSh^2))
Hotelling_correction=(2*(n-1))/(n-2)
ell_OLSh=ellipse(cov_matrix_OLSh,centre=c(slope_OLSh,intercept_OLSh),npoints=1000, t = sqrt(qf(conf.level,2,n-2)*Hotelling_correction))
F1=matrix(nrow=1,ncol=2,c(slope_OLSh-1,intercept_OLSh))
F2=solve(cov_matrix_OLSh)
F3=t(F1)
Fell=F1%*%F2%*%F3
alpha=1-conf.level
CI_slope_OLSh_1_xy=slope_OLSh-qt(1-alpha/2,n-2)*S_slope_OLSh
CI_slope_OLSh_2_xy=slope_OLSh+qt(1-alpha/2,n-2)*S_slope_OLSh
CI_intercept_OLSh_1_xy=intercept_OLSh-qt(1-alpha/2,n-2)*S_intercept_OLSh
CI_intercept_OLSh_2_xy=intercept_OLSh+qt(1-alpha/2,n-2)*S_intercept_OLSh
pval_slope_OLSh=(1-pt(abs(slope_OLSh-1)/S_slope_OLSh,n-2))*2
pval_intercept_OLSh=(1-pt(abs(intercept_OLSh)/S_intercept_OLSh,n-2))*2
pval_ell_OLSh=1-pf(F1%*%F2%*%F3/Hotelling_correction,2,n-2)

ell_OLSh[,2]=-ell_OLSh[,2]/slope_OLSh
ell_OLSh[,1]=1/ell_OLSh[,1]
intercept_OLSh=-intercept_OLSh/slope_OLSh
CI_intercept_OLSh_1=-CI_intercept_OLSh_2_xy/slope_OLSh
CI_intercept_OLSh_2=-CI_intercept_OLSh_1_xy/slope_OLSh
slope_OLSh=1/slope_OLSh
CI_slope_OLSh_1=1/CI_slope_OLSh_2_xy
CI_slope_OLSh_2=1/CI_slope_OLSh_1_xy

rownames=c("Intercept","Slope","Joint")
name1=paste("Lower ",conf.level*100,"%CI",sep="")
name2=paste("Upper ",conf.level*100,"%CI",sep="")
colnames=c("H0","Estimate",name1,name2,"pvalue")
Table.OLSh=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))

Table.OLSh$H0=c("0","1","(0,1)")
Table.OLSh$Estimate=c(intercept_OLSh,slope_OLSh,"")
Table.OLSh[,name1]=c(CI_intercept_OLSh_1,CI_slope_OLSh_1,"")
Table.OLSh[,name2]=c(CI_intercept_OLSh_2,CI_slope_OLSh_2,"")
Table.OLSh$pvalue=c(pval_intercept_OLSh,pval_slope_OLSh,pval_ell_OLSh)

results=list(ell_OLSh,Table.OLSh)
names(results)=c("Ellipse.OLSh","Estimate.OLSh")
return(results)
}
