#' @export
BLS.fit <-
function(data=NULL,xcol=1,ycol=2,var.x=NULL,var.y=NULL,ratio.var=NULL,conf.level=0.95)
{

if(length(conf.level) > 1) stop("Only one confidence level is allowed.")
if(conf.level <= 0 | conf.level >= 1) stop("The confidence level must be between 0 and 1 (excluded)")
if(length(ratio.var) > 1) stop("Only one value is allowed in ratio.var")
if(!is.null(ratio.var)) {if(ratio.var <= 0 | !is.numeric(ratio.var)) stop("ratio.var must be numeric and positive")}
if(length(var.x) > 1) stop("var.x must be a number (greater than 0)")
if(!is.null(var.x)) {if(var.x <= 0 | !is.numeric(var.x)) stop("var.x must be a number (greater than 0)")}
if(length(var.y) > 1) stop("var.y must be a number (greater than 0)")
if(!is.null(var.y)) {if(var.y <= 0 | !is.numeric(var.y)) stop("var.y must be a number (greater than 0)")}


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

alpha.conf=1-conf.level
CI_slope_BLS_1=slope_BLS-qt(1-alpha.conf/2,n-2)*S_slope_BLS
CI_slope_BLS_2=slope_BLS+qt(1-alpha.conf/2,n-2)*S_slope_BLS
CI_intercept_BLS_1=intercept_BLS-qt(1-alpha.conf/2,n-2)*S_intercept_BLS
CI_intercept_BLS_2=intercept_BLS+qt(1-alpha.conf/2,n-2)*S_intercept_BLS
pval_slope_BLS=(1-pt(abs(slope_BLS-1)/S_slope_BLS,n-2))*2
pval_intercept_BLS=(1-pt(abs(intercept_BLS)/S_intercept_BLS,n-2))*2

rownames=c("Intercept","Slope")
name_CI1=paste("Lower ",conf.level*100,"%CI",sep="")
name_CI2=paste("Upper ",conf.level*100,"%CI",sep="")
colnames=c("H0","Estimate","Std Error",name_CI1,name_CI2,"pvalue")
Table.BLS=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))

Table.BLS$H0=c(0,1)
Table.BLS$Estimate=c(intercept_BLS,slope_BLS)
Table.BLS$'Std Error'=c(S_intercept_BLS,S_slope_BLS)
Table.BLS[,name_CI1]=c(CI_intercept_BLS_1,CI_slope_BLS_1)
Table.BLS[,name_CI2]=c(CI_intercept_BLS_2,CI_slope_BLS_2)
Table.BLS$pvalue=c(pval_intercept_BLS,pval_slope_BLS)

return(Table.BLS)
}
