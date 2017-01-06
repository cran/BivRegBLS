#' @export
DR <-
function(data=NULL,xcol=1,ycol=2,ratio.var=NULL,conf.level=0.95)
{

if(length(conf.level) > 1) stop("Only one confidence level is allowed.")
if(conf.level <= 0 | conf.level >= 1) stop("The confidence level must be between 0 and 1 (excluded)")
if(length(ratio.var) > 1) stop("Only one value is allowed in ratio.var")
if(!is.null(ratio.var)) {if(ratio.var <= 0 | !is.numeric(ratio.var)) stop("ratio.var must be numeric and positive")}

res=desc.stat(data=data,xcol=xcol,ycol=ycol)
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
varx=res$statistics$'Var X'
vary=res$statistics$'Var Y'

if(is.null(ratio.var))
{
	if (nx==1 | ny==1)
	{
	text1="Lambda (ratio of measurement errors variances) is not given and not estimatable (unreplicated data) \n"
	text2="It is therefore not possible to estimate a regression line \n"
	text3="Try the functions 'FullCIs' and 'GraphFullCIs' to estimate all the regression lines from OLSv to OLSh \n"
	textall=paste(text1,text2,text3)
	cat(textall)
	return()
	}
	res=lambdas(data=data,xcol=xcol,ycol=ycol)
	lambda=res$lambdas.Table$Estimate[2]
} else {
lambda=ratio.var
varx=NA
vary=NA
}

slope_DR=(Syy-lambda*Sxx+sqrt((Syy-lambda*Sxx)^2+4*lambda*Sxy^2))/(2*Sxy)
intercept_DR=ymean-xmean*slope_DR
tetas=atan(slope_DR/sqrt(lambda))
alpha=1-conf.level
arg_asin=2*qt(1-alpha/2,n-2)/sqrt(n-2)*sqrt((lambda*(Syy*Sxx-Sxy^2))/((Syy-lambda*Sxx)^2+4*lambda*Sxy^2))
if(arg_asin >= -1 | arg_asin <= 1) phi=asin(arg_asin)/2
if(arg_asin >= -1 | arg_asin <= 1) CI_slope_DR_exact_1=sqrt(lambda)*tan(tetas-phi) else CI_slope_DR_exact_1=NA
if(arg_asin >= -1 | arg_asin <= 1) CI_slope_DR_exact_2=sqrt(lambda)*tan(tetas+phi) else CI_slope_DR_exact_2=NA

S_slope_DR=sqrt((Sxx*Syy-Sxy^2)/(n*(Sxy/slope_DR)^2))
S_intercept_DR=sqrt(xmean^2*S_slope_DR^2+(slope_DR^2*varx/nx+vary/ny)/n)
CI_slope_DR_approx_1=slope_DR-qt(1-alpha/2,n-2)*S_slope_DR
CI_slope_DR_approx_2=slope_DR+qt(1-alpha/2,n-2)*S_slope_DR
CI_intercept_DR_approx_1=intercept_DR-qt(1-alpha/2,n-2)*S_intercept_DR
CI_intercept_DR_approx_2=intercept_DR+qt(1-alpha/2,n-2)*S_intercept_DR
pval_slope_DR_approx=(1-pt(abs(slope_DR-1)/S_slope_DR,n-2))*2
pval_intercept_DR_approx=(1-pt(abs(intercept_DR)/S_intercept_DR,n-2))*2

cov_slope_intercept_DR=-S_slope_DR^2*xmean
cov_matrix_DR=matrix(nrow=2,ncol=2,c(S_slope_DR^2,cov_slope_intercept_DR,cov_slope_intercept_DR,S_intercept_DR^2))
Hotelling_correction=(2*(n-1))/(n-2)
if(sum(is.na(cov_matrix_DR))>0) ell_DR=NULL else ell_DR=ellipse(cov_matrix_DR,centre=c(slope_DR,intercept_DR),npoints=1000, t = sqrt(qf(conf.level,2,n-2)*Hotelling_correction))

F1=matrix(nrow=1,ncol=2,c(slope_DR-1,intercept_DR))
F2=solve(cov_matrix_DR)
F3=t(F1)
Fell=F1%*%F2%*%F3
pval_ell_DR=1-pf(F1%*%F2%*%F3/Hotelling_correction,2,n-2)

rownames=c("Intercept (approximate)","Slope (exact)","Slope (approximate)","Joint")
name_CI1=paste("Lower ",conf.level*100,"%CI",sep="")
name_CI2=paste("Upper ",conf.level*100,"%CI",sep="")
colnames=c("H0","Estimate",name_CI1,name_CI2,"pvalue")
Table.DR=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))

Table.DR$H0=c(0,1,1,"(0,1)")
Table.DR$Estimate=c(intercept_DR,slope_DR,slope_DR,"")
Table.DR[,name_CI1]=c(CI_intercept_DR_approx_1,CI_slope_DR_exact_1,CI_slope_DR_approx_1,"")
Table.DR[,name_CI2]=c(CI_intercept_DR_approx_2,CI_slope_DR_exact_2,CI_slope_DR_approx_2,"")
Table.DR$pvalue=c(pval_intercept_DR_approx,"",pval_slope_DR_approx,pval_ell_DR)

if(sum(is.na(cov_matrix_DR))>0) Table.DR=Table.DR[which(rownames(Table.DR) != "Joint"),]

results=list(ell_DR,Table.DR)
names(results)=c("Ellipse.DR","Estimate.DR")
return(results)
}
