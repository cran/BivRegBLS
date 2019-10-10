#' @export
CBLS.fit <-
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
x=res$Xi
y=res$Yi
Xij=res$Xij
Yik=res$Yik
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

d=y-x
m=(x+y)/2
y=d
x=m
rhoMD=(lambda-1)/(lambda+1)
lambdaxy=4
xmean=mean(x)
xcentre=x-xmean
Sxx=sum((x-xmean)^2)
ymean=mean(y)
ycentre=y-ymean
Syy=sum((y-ymean)^2)
Sxy=sum(xcentre*ycentre)

slope_CBLS=(Syy-4*Sxx+sqrt((Syy-4*Sxx)^2-4*(Sxy-2*rhoMD*Sxx)*(2*rhoMD*Syy-4*Sxy)))/(2*(Sxy-2*rhoMD*Sxx))
intercept_CBLS=ymean-xmean*slope_CBLS

varXX=(varX+varY)/4
varYY=varX+varY
varX=varXX
varY=varYY
cov_CBLS=rhoMD*sqrt(varX)*sqrt(varY)

R=matrix(nrow=2,ncol=2)
g=matrix(nrow=2,ncol=1)
#diff=matrix(nrow=2,ncol=1)
#repeat {
	ei=y-(x*slope_CBLS+intercept_CBLS)
#	ei2=ei^2
	Sei2=varY+slope_CBLS^2*varX-2*slope_CBLS*cov_CBLS
	R1i=1/Sei2
	R2i=x*R1i
	R4i=x^2*R1i	
	R[1,1]=sum(R1i)
	R[1,2]=sum(R2i)
	R[2,1]=sum(R2i)
	R[2,2]=sum(R4i)
	g1i=y*R1i
	g2i=x*y*R1i+(ei*R1i)^2*(slope_CBLS*varX-cov_CBLS)
	g[1,1]=sum(g1i)
	g[2,1]=sum(g2i)
	b=solve(R)%*%g
#	diff=b0-b
#	b0=b0-diff/2
#	if ((abs(diff[1,1])<0.000001) & (abs(diff[2,1])<0.000001))
#	break()
#	}
slope_CBLS=b[2,1]
intercept_CBLS=b[1,1]
s2i=sum(((y-(x*slope_CBLS+intercept_CBLS))^2)/Sei2)/(n-2)
varcovari=solve(R)*s2i
S_slope_CBLS=sqrt(varcovari[2,2])
S_intercept_CBLS=sqrt(varcovari[1,1])

alpha.conf=1-conf.level
CI_slope_CBLS_1=slope_CBLS-qt(1-alpha.conf/2,n-2)*S_slope_CBLS
CI_slope_CBLS_2=slope_CBLS+qt(1-alpha.conf/2,n-2)*S_slope_CBLS
CI_intercept_CBLS_1=intercept_CBLS-qt(1-alpha.conf/2,n-2)*S_intercept_CBLS
CI_intercept_CBLS_2=intercept_CBLS+qt(1-alpha.conf/2,n-2)*S_intercept_CBLS
pval_slope_CBLS=(1-pt(abs(slope_CBLS-0)/S_slope_CBLS,n-2))*2
pval_intercept_CBLS=(1-pt(abs(intercept_CBLS)/S_intercept_CBLS,n-2))*2

rownames=c("Intercept","Slope")
name_CI1=paste("Lower ",conf.level*100,"%CI",sep="")
name_CI2=paste("Upper ",conf.level*100,"%CI",sep="")
colnames=c("H0","Estimate","Std Error",name_CI1,name_CI2,"pvalue")
Table.CBLS=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))

Table.CBLS$H0=c("0","0")
Table.CBLS$Estimate=c(intercept_CBLS,slope_CBLS)
Table.CBLS$'Std Error'=c(S_intercept_CBLS,S_slope_CBLS)
Table.CBLS[,name_CI1]=c(CI_intercept_CBLS_1,CI_slope_CBLS_1)
Table.CBLS[,name_CI2]=c(CI_intercept_CBLS_2,CI_slope_CBLS_2)
Table.CBLS$pvalue=c(pval_intercept_CBLS,pval_slope_CBLS)

return(Table.CBLS)
}
