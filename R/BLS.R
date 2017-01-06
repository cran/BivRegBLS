#' @export
BLS <-
function(data=NULL,xcol=1,ycol=2,var.x=NULL,var.y=NULL,df.var.x=Inf,df.var.y=Inf,ratio.var=NULL,conf.level=0.95,pred.level=0.95,npoints=1000,qx=1,qy=1,xpred=NULL)
{

if(length(conf.level) > 1) stop("Only one confidence level is allowed.")
if(length(pred.level) > 1) stop("Only one predictive level is allowed.")
if(conf.level <= 0 | conf.level >= 1) stop("The confidence level must be between 0 and 1 (excluded)")
if(pred.level <= 0 | pred.level >= 1) stop("The predictive level must be between 0 and 1 (excluded)")
if(length(ratio.var) > 1) stop("Only one value is allowed in ratio.var")
if(!is.null(ratio.var)) {if(ratio.var <= 0 | !is.numeric(ratio.var)) stop("ratio.var must be numeric and positive")}
if(length(var.x) > 1) stop("var.x must be a number (greater than 0)")
if(!is.null(var.x)) {if(var.x <= 0 | !is.numeric(var.x)) stop("var.x must be a number (greater than 0)")}
if(length(var.y) > 1) stop("var.y must be a number (greater than 0)")
if(!is.null(var.y)) {if(var.y <= 0 | !is.numeric(var.y)) stop("var.y must be a number (greater than 0)")}
if(length(df.var.x) > 1) stop("df.var.x must be a number (greater than 0)")
if(df.var.x <= 0 | !is.numeric(df.var.x)) stop("df.var.x must be a number (greater than 0)")
if(length(df.var.y) > 1) stop("df.var.y must be a number (greater than 0)")
if(df.var.y <= 0 | !is.numeric(df.var.y)) stop("df.var.y must be a number (greater than 0)")
if(length(npoints) > 1) stop("npoints must be an integer (at least 10)")
if(floor(npoints) != ceiling(npoints) | npoints < 10 | !is.numeric(npoints)) stop("npoints must be an integer (at least 10)")
if(length(qx) > 1) stop("qx must be an integer (at least 1)")
if(floor(qx) != ceiling(qx) | qx < 1 | !is.numeric(qx)) stop("qx must be an integer (at least 1)")
if(length(qy) > 1) stop("qy must be an integer (at least 1)")
if(floor(qy) != ceiling(qy) | qy < 1 | !is.numeric(qy)) stop("qy must be an integer (at least 1)")
if(!is.null(xpred) & !is.numeric(xpred)) stop("xpred must be numeric")

if(is.null(var.x) & df.var.x != Inf) warning("df.var.x is ignored as var.x is not given")
if(is.null(var.y) & df.var.y != Inf) warning("df.var.y is ignored as var.y is not given")

res=desc.stat(data=data,xcol=xcol,ycol=ycol)
Xij=res$Xij
Yik=res$Yik
x=res$Xi
y=res$Yi
xmean=res$statistics$Xmean
ymean=res$statistics$Ymean
Sxx=res$statistics$Sxx
Syy=res$statistics$Syy
Sxy=res$statistics$Sxy
n=res$statistics$N
nx=res$statistics$nx
ny=res$statistics$ny
Nx=res$statistics$df.var.x
Ny=res$statistics$df.var.y
nxi=res$nxi
nyi=res$nyi
if(is.null(var.x)) variances_x=res$variances_x else variances_x=var.x
if(is.null(var.y)) variances_y=res$variances_y else variances_y=var.y

if((!is.null(var.x) & !is.null(var.y) & !is.null(ratio.var)) | ((nx>1 | ny>1) & !is.null(ratio.var))) warning("The value of ratio.var is ignored if var.x and var.y are given, or if replicates are available.")

if(nx > 1 & !is.null(var.x)) warning("The measurement error variance on X is given by var.x, it will not be calculated from the X replicates.")
if(ny > 1 & !is.null(var.y)) warning("The measurement error variance on Y is given by var.y, it will not be calculated from the Y replicates.")

if(nx==1 & ny==1)
{
	if((is.null(var.x) | is.null(var.y)) & is.null(ratio.var))
	{
	text1="At least one measurement errors variance is not given and not estimatable (unreplicated data) \n"
	text2="Lambda (ratio of measurement error variances) is not estimatable \n"
	text3="It is therefore not possible to estimate a regression line \n"
	text4="Try the functions 'FullCIs' and 'GraphFullCIs' to estimate all the regression lines from OLSv to OLSh \n"
	textall=paste(text1,text2,text3,text4)
	cat(textall)
	return()
	}
	if ((is.null(var.x) + is.null(var.y) >= 1) & !is.null(ratio.var))
	{
		if((is.null(var.x) & !is.null(var.y))) warning("The value of var.y is ignored as lambda (ratio.var) is given and var.x is not estimatable (unreplicated data)")
		if((!is.null(var.x) & is.null(var.y))) warning("The value of var.x is ignored as lambda (ratio.var) is given and var.y is not estimatable (unreplicated data)")
	lambda=ratio.var
	varX=rep(1,n)
	varY=rep(ratio.var,n)
	varx=NA
	vary=NA
	}
	if(!is.null(var.x) & !is.null(var.y))
	{
	lambda=var.y/var.x
	varX=rep(var.x,n)
	varY=rep(var.y,n)
	varx=var.x
	vary=var.y
	df.var.x=df.var.x
	df.var.y=df.var.y
	}
}

if(nx==1 & ny>1)
{
	if(is.null(var.x) & is.null(ratio.var))
	{
	text1="The measurement errors variance on X is not given and not estimatable (unreplicated data) \n"
	text2="Lambda (ratio of measurement error variances) is not estimatable \n"
	text3="It is therefore not possible to estimate a regression line \n"
	text4="Try the functions 'FullCIs' and 'GraphFullCIs' to estimate all the regression lines from OLSv to OLSh \n"
	textall=paste(text1,text2,text3,text4)
	cat(textall)
	return()
	}
	if(is.null(var.x) & !is.null(ratio.var))
	{
	lambda=ratio.var
	varX=rep(1,n)
	varY=rep(ratio.var,n)
	varx=NA
	vary=NA
	}
	if(!is.null(var.x))
	{
	res=suppressWarnings(desc.stat(data=data,xcol=xcol,ycol=ycol))
	if(is.null(var.y)) vary=res$statistics$'Var Y' else vary=var.y
	varx=var.x
	lambda=(vary/ny)/varx
	varX=rep(varx,n)
	varY=rep(vary/ny,n)
	df.var.x=df.var.x
	if(is.null(var.y)) df.var.y=Ny else df.var.y=df.var.y
	}
}

if(nx>1 & ny==1)
{
	if(is.null(var.y) & is.null(ratio.var))
	{
	text1="The measurement errors variance on Y is not given and not estimatable (unreplicated data) \n"
	text2="Lambda (ratio of measurement error variances) is not estimatable \n"
	text3="It is therefore not possible to estimate a regression line \n"
	text4="Try the functions 'FullCIs' and 'GraphFullCIs' to estimate all the regression lines from OLSv to OLSh \n"
	textall=paste(text1,text2,text3,text4)
	cat(textall)
	return()
	}
	if(is.null(var.y) & !is.null(ratio.var))
	{
	lambda=ratio.var
	varX=rep(1,n)
	varY=rep(ratio.var,n)
	varx=NA
	vary=NA
	}
	if(!is.null(var.y))
	{
	res=suppressWarnings(desc.stat(data=data,xcol=xcol,ycol=ycol))
	if(is.null(var.x)) varx=res$statistics$'Var X' else varx=var.x
	vary=var.y
	lambda=vary/(varx/nx)
	varX=rep(varx/nx,n)
	varY=rep(vary,n)
	if(is.null(var.x)) df.var.x=Nx else df.var.x=df.var.x
	df.var.y=df.var.y
	}
}

if(nx>1 & ny>1)
{
	res=suppressWarnings(desc.stat(data=data,xcol=xcol,ycol=ycol))
	if(is.null(var.x)) varx=res$statistics$'Var X' else varx=var.x
	if(is.null(var.y)) vary=res$statistics$'Var Y' else vary=var.y
	lambda=(vary/ny)/(varx/nx)
	varX=rep(varx/nx,n)
	varY=rep(vary/ny,n)
	if(is.null(var.x)) df.var.x=Nx else df.var.x=df.var.x
	if(is.null(var.y)) df.var.y=Ny else df.var.y=df.var.y
}


slope_BLS=(Syy-lambda*Sxx+sqrt((Syy-lambda*Sxx)^2+4*lambda*Sxy^2))/(2*Sxy)
intercept_BLS=ymean-xmean*slope_BLS

R=matrix(nrow=2,ncol=2)
g=matrix(nrow=2,ncol=1)
#diff=matrix(nrow=2,ncol=1)
#repeat {
	ei=y-(x*slope_BLS+intercept_BLS)
#	ei2=ei^2
	Sei2=varY+slope_BLS^2*varX
	R1i=1/Sei2
	R2i=x*R1i
	R4i=x^2*R1i	
	R[1,1]=sum(R1i)
	R[1,2]=sum(R2i)
	R[2,1]=sum(R2i)
	R[2,2]=sum(R4i)
	g1i=y*R1i
	g2i=x*y*R1i+slope_BLS*(ei*R1i)^2*varX
	g[1,1]=sum(g1i)
	g[2,1]=sum(g2i)
	b=solve(R)%*%g
#	diff=b0-b
#	b0=b0-diff/2
#	if ((abs(diff[1,1])<0.000001) & (abs(diff[2,1])<0.000001))
#	break()
#	}
slope_BLS=b[2,1]
intercept_BLS=b[1,1]
s2i=sum(((y-(x*slope_BLS+intercept_BLS))^2)/Sei2)/(n-2)
varcovari=solve(R)*s2i
S_slope_BLS=sqrt(varcovari[2,2])
S_intercept_BLS=sqrt(varcovari[1,1])
cov_slope_intercept_BLS=varcovari[1,2]
cov_matrix_BLS=matrix(nrow=2,ncol=2,c(S_slope_BLS^2,cov_slope_intercept_BLS,cov_slope_intercept_BLS,S_intercept_BLS^2))
Hotelling_correction=(2*(n-1))/(n-2)
ell_BLS=ellipse(cov_matrix_BLS,centre=c(slope_BLS,intercept_BLS),npoints=npoints, t = sqrt(qf(conf.level,2,n-2)*Hotelling_correction))

alpha.conf=1-conf.level
CI_slope_BLS_1=slope_BLS-qt(1-alpha.conf/2,n-2)*S_slope_BLS
CI_slope_BLS_2=slope_BLS+qt(1-alpha.conf/2,n-2)*S_slope_BLS
CI_intercept_BLS_1=intercept_BLS-qt(1-alpha.conf/2,n-2)*S_intercept_BLS
CI_intercept_BLS_2=intercept_BLS+qt(1-alpha.conf/2,n-2)*S_intercept_BLS
pval_slope_BLS=(1-pt(abs(slope_BLS-1)/S_slope_BLS,n-2))*2
pval_intercept_BLS=(1-pt(abs(intercept_BLS)/S_intercept_BLS,n-2))*2

F1=matrix(nrow=1,ncol=2,c(slope_BLS-1,intercept_BLS))
F2=solve(cov_matrix_BLS)
F3=t(F1)
Fell=F1%*%F2%*%F3
pval_ell_BLS=1-pf(F1%*%F2%*%F3/Hotelling_correction,2,n-2)

rownames=c("Intercept","Slope","Joint")
name_CI1=paste("Lower ",conf.level*100,"%CI",sep="")
name_CI2=paste("Upper ",conf.level*100,"%CI",sep="")
colnames=c("H0","Estimate","Std Error",name_CI1,name_CI2,"pvalue")
Table.BLS=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))

Table.BLS$H0=c("0","1","(0,1)")
Table.BLS$Estimate=c(intercept_BLS,slope_BLS,"")
Table.BLS$'Std Error'=c(S_intercept_BLS,S_slope_BLS,"")
Table.BLS[,name_CI1]=c(CI_intercept_BLS_1,CI_slope_BLS_1,"")
Table.BLS[,name_CI2]=c(CI_intercept_BLS_2,CI_slope_BLS_2,"")
Table.BLS$pvalue=c(pval_intercept_BLS,pval_slope_BLS,pval_ell_BLS)

xx=c(seq(min(x),max(x),length.out=npoints),xpred)
xx.pred=xx*slope_BLS+intercept_BLS
var.line=S_intercept_BLS^2+S_slope_BLS^2*xx^2+2*xx*cov_slope_intercept_BLS

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

W.BLS=varY[1]+varX[1]*slope_BLS^2
alpha.pred=1-pred.level

data.pred[,CI1]=xx.pred-qt(1-alpha.conf/2,n-2)*sqrt(var.line)
data.pred[,CI2]=xx.pred+qt(1-alpha.conf/2,n-2)*sqrt(var.line)

var_y_single=vary
var_x_single=varx
var_PI=var.line+var_y_single+var_x_single*slope_BLS^2
Numerator=(var.line+var_y_single+slope_BLS^2*var_x_single)^2
Denominator=var.line^2/(n-2)+var_y_single^2/df.var.y+(slope_BLS^2*var_x_single)^2/df.var.x
dfPI=Numerator/Denominator
data.pred[,PI1]=xx.pred-qt(1-alpha.pred/2,dfPI)*sqrt(var_PI)
data.pred[,PI2]=xx.pred+qt(1-alpha.pred/2,dfPI)*sqrt(var_PI)

var_y_q=vary/qy
var_x_q=varx/qx
var_GI=var.line+var_y_q+var_x_q*slope_BLS^2
Numerator=(var.line+var_y_q+slope_BLS^2*var_x_q)^2
Denominator=var.line^2/(n-2)+var_y_q^2/df.var.y+(slope_BLS^2*var_x_q)^2/df.var.x
dfGI=Numerator/Denominator
data.pred[,GI1]=xx.pred-qt(1-alpha.pred/2,dfGI)*sqrt(var_GI)
data.pred[,GI2]=xx.pred+qt(1-alpha.pred/2,dfGI)*sqrt(var_GI)

data.pred[,CB1]=xx.pred-sqrt(qf(conf.level,2,n-2)*Hotelling_correction)*sqrt((W.BLS/(n*sum(x^2)-sum(x)^2)*(sum(x^2)-2*xx*sum(x)+n*xx^2))*s2i)
data.pred[,CB2]=xx.pred+sqrt(qf(conf.level,2,n-2)*Hotelling_correction)*sqrt((W.BLS/(n*sum(x^2)-sum(x)^2)*(sum(x^2)-2*xx*sum(x)+n*xx^2))*s2i)

results=list(Xij,Yik,x,y,nxi,nyi,variances_x,variances_y,lambda,ell_BLS,Table.BLS,data.pred[1:npoints,],data.pred[-(1:npoints),])
names(results)=c("Xij","Yik","Xi","Yi","nxi","nyi","variances_x","variances_y","Lambda.XY","Ellipse.BLS","Estimate.BLS","Pred.BLS","xpred.BLS")
class(results)="BLS"
return(results)
}
