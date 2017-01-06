#' @importFrom ellipse ellipse
#' @export
OLSv <-
function(data=NULL,xcol=1,ycol=2,conf.level=0.95,pred.level=0.95,npoints=1000,q=1,xpred=NULL)
{

if(length(conf.level) > 1) stop("Only one confidence level is allowed.")
if(length(pred.level) > 1) stop("Only one predictive level is allowed.")
if(conf.level <= 0 | conf.level >= 1) stop("The confidence level must be between 0 and 1 (excluded)")
if(pred.level <= 0 | pred.level >= 1) stop("The predictive level must be between 0 and 1 (excluded)")
if(length(npoints) > 1) stop("npoints must be an integer (at least 10)")
if(floor(npoints) != ceiling(npoints) | npoints < 10 | !is.numeric(npoints)) stop("npoints must be an integer (at least 10)")
if(length(q) > 1) stop("q must be an integer (at least 1)")
if(floor(q) != ceiling(q) | q < 1 | !is.numeric(q)) stop("q must be an integer (at least 1)")
if(!is.null(xpred) & !is.numeric(xpred)) stop("xpred must be numeric")

res=desc.stat(data=data,xcol=xcol,ycol=ycol)
x=res$Xi
y=res$Yi
xmean=res$statistics$Xmean
ymean=res$statistics$Ymean
Sxx=res$statistics$Sxx
Syy=res$statistics$Syy
Sxy=res$statistics$Sxy
n=res$statistics$N

slope_OLSv=Sxy/Sxx
intercept_OLSv=ymean-xmean*slope_OLSv
S2_OLSv=sum((y-intercept_OLSv-slope_OLSv*x)^2)/(n-2)
S_slope_OLSv=sqrt(S2_OLSv/Sxx)
S_intercept_OLSv=sqrt(S2_OLSv*(1/n+xmean^2/Sxx))
cov_slope_intercept_OLSv=-S2_OLSv*xmean/Sxx
cov_matrix_OLSv=matrix(nrow=2,ncol=2,c(S_slope_OLSv^2,cov_slope_intercept_OLSv,cov_slope_intercept_OLSv,S_intercept_OLSv^2))
Hotelling_correction=(2*(n-1))/(n-2)
ell_OLSv=ellipse(cov_matrix_OLSv,centre=c(slope_OLSv,intercept_OLSv),npoints=npoints, t = sqrt(qf(conf.level,2,n-2)*Hotelling_correction))
F1=matrix(nrow=1,ncol=2,c(slope_OLSv-1,intercept_OLSv))
F2=solve(cov_matrix_OLSv)
F3=t(F1)
Fell=F1%*%F2%*%F3
alpha=1-conf.level
CI_slope_OLSv_1=slope_OLSv-qt(1-alpha/2,n-2)*S_slope_OLSv
CI_slope_OLSv_2=slope_OLSv+qt(1-alpha/2,n-2)*S_slope_OLSv
CI_intercept_OLSv_1=intercept_OLSv-qt(1-alpha/2,n-2)*S_intercept_OLSv
CI_intercept_OLSv_2=intercept_OLSv+qt(1-alpha/2,n-2)*S_intercept_OLSv
pval_slope_OLSv=(1-pt(abs(slope_OLSv-1)/S_slope_OLSv,n-2))*2
pval_intercept_OLSv=(1-pt(abs(intercept_OLSv)/S_intercept_OLSv,n-2))*2
pval_ell_OLSv=1-pf(F1%*%F2%*%F3/Hotelling_correction,2,n-2)

rownames=c("Intercept","Slope","Joint")
name1=paste("Lower ",conf.level*100,"%CI",sep="")
name2=paste("Upper ",conf.level*100,"%CI",sep="")
colnames=c("H0","Estimate","Std Error",name1,name2,"pvalue")
Table.OLSv=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))

Table.OLSv$H0=c("0","1","(0,1)")
Table.OLSv$Estimate=c(intercept_OLSv,slope_OLSv,"")
Table.OLSv$'Std Error'=c(S_intercept_OLSv,S_slope_OLSv,"")
Table.OLSv[,name1]=c(CI_intercept_OLSv_1,CI_slope_OLSv_1,"")
Table.OLSv[,name2]=c(CI_intercept_OLSv_2,CI_slope_OLSv_2,"")
Table.OLSv$pvalue=c(pval_intercept_OLSv,pval_slope_OLSv,pval_ell_OLSv)

xx=c(seq(min(x),max(x),length.out=npoints),xpred)
xx.pred=xx*slope_OLSv+intercept_OLSv

CI1=paste(conf.level*100,"% CI Lower",sep="")
CI2=paste(conf.level*100,"% CI Upper",sep="")
CB1=paste(conf.level*100,"% CB Lower",sep="")
CB2=paste(conf.level*100,"% CB Upper",sep="")
PI1=paste(pred.level*100,"% PI Lower",sep="")
PI2=paste(pred.level*100,"% PI Upper",sep="")
GI1=paste(pred.level*100,"% GI Lower",sep="")
GI2=paste(pred.level*100,"% GI Upper",sep="")

colnames=c("X0","Ypred",CI1,CI2,PI1,PI2,GI1,GI2,CB1,CB2)
if(is.null(xpred)) rownames=1:npoints else rownames=c(1:npoints,paste("xpred",1:length(xpred)))
data.pred=as.data.frame(matrix(nrow=npoints+length(xpred),ncol=length(colnames),dimnames=list(rownames,colnames)))
data.pred$X0=xx
data.pred$Ypred=xx.pred
alpha.pred=1-pred.level

data.pred[,CI1]=xx.pred-qt(1-alpha/2,n-2)*sqrt(S2_OLSv)*sqrt(1/n+(xx-xmean)^2/Sxx)
data.pred[,CI2]=xx.pred+qt(1-alpha/2,n-2)*sqrt(S2_OLSv)*sqrt(1/n+(xx-xmean)^2/Sxx)

data.pred[,PI1]=xx.pred-qt(1-alpha.pred/2,n-2)*sqrt(S2_OLSv)*sqrt(1+1/n+(xx-xmean)^2/Sxx)
data.pred[,PI2]=xx.pred+qt(1-alpha.pred/2,n-2)*sqrt(S2_OLSv)*sqrt(1+1/n+(xx-xmean)^2/Sxx)

data.pred[,GI1]=xx.pred-qt(1-alpha.pred/2,n-2)*sqrt(S2_OLSv)*sqrt(1/q+1/n+(xx-xmean)^2/Sxx)
data.pred[,GI2]=xx.pred+qt(1-alpha.pred/2,n-2)*sqrt(S2_OLSv)*sqrt(1/q+1/n+(xx-xmean)^2/Sxx)

data.pred[,CB1]=xx.pred-sqrt(qf(conf.level,2,n-2)*Hotelling_correction)*sqrt(S2_OLSv)*sqrt(1/n+(xx-xmean)^2/Sxx)
data.pred[,CB2]=xx.pred+sqrt(qf(conf.level,2,n-2)*Hotelling_correction)*sqrt(S2_OLSv)*sqrt(1/n+(xx-xmean)^2/Sxx)

results=list(ell_OLSv,Table.OLSv,data.pred[1:npoints,],data.pred[-(1:npoints),])
names(results)=c("Ellipse.OLSv","Estimate.OLSv","Pred.OLSv","xpred.OLSv")
return(results)
}
