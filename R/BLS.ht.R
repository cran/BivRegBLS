#' @export
BLS.ht <-
function(data=NULL,xcol=NULL,ycol=NULL,var.x.col=NULL,var.y.col=NULL,var.x.formula=NULL,var.y.formula=NULL,nxi.col=NULL,nyi.col=NULL,df.var.x.col=NULL,df.var.y.col=NULL,conf.level=0.95,pred.level=0.95,npoints=1000,qx=1,qy=1,xpred=NULL,tol=1e-05)
{

if(length(conf.level) > 1) stop("Only one confidence level is allowed.")
if(length(pred.level) > 1) stop("Only one predictive level is allowed.")
if(conf.level <= 0 | conf.level >= 1) stop("The confidence level must be between 0 and 1 (excluded)")
if(pred.level <= 0 | pred.level >= 1) stop("The predictive level must be between 0 and 1 (excluded)")
if(length(npoints) > 1) stop("npoints must be an integer (at least 10)")
if(floor(npoints) != ceiling(npoints) | npoints < 10 | !is.numeric(npoints)) stop("npoints must be an integer (at least 10)")
if(length(qx) > 1) stop("qx must be an integer (at least 1)")
if(floor(qx) != ceiling(qx) | qx < 1 | !is.numeric(qx)) stop("qx must be an integer (at least 1)")
if(length(qy) > 1) stop("qy must be an integer (at least 1)")
if(floor(qy) != ceiling(qy) | qy < 1 | !is.numeric(qy)) stop("qy must be an integer (at least 1)")
if(!is.null(xpred) & !is.numeric(xpred)) stop("xpred must be numeric")
if(length(var.x.formula) > 1) stop("var.x.formula must be character (length = 1)")
if(!is.null(var.x.formula) & !is.character(var.x.formula)) stop("var.x.formula must be character")
if(length(var.y.formula) > 1) stop("var.y.formula must be character (length = 1)")
if(!is.null(var.y.formula) & !is.character(var.y.formula)) stop("var.y.formula must be character")
if(length(tol) > 1 | !is.numeric(tol)) stop("tol must be one numeric value")
if(tol <= 0) stop("tol must be one numeric value greater than 0")

if(is.null(colnames(data))) colnames(data)=paste0("V",1:ncol(data))
if(!is.null(var.x.col) & is.numeric(var.x.col)) var.x.col=colnames(data)[var.x.col]
if(!is.null(var.y.col) & is.numeric(var.y.col)) var.y.col=colnames(data)[var.y.col]
if(!is.null(df.var.x.col) & is.numeric(df.var.x.col)) df.var.x.col=colnames(data)[df.var.x.col]
if(!is.null(df.var.y.col) & is.numeric(df.var.y.col)) df.var.y.col=colnames(data)[df.var.y.col]
if(!is.null(nxi.col) & is.numeric(nxi.col)) nxi.col=colnames(data)[nxi.col]
if(!is.null(nyi.col) & is.numeric(nyi.col)) nyi.col=colnames(data)[nyi.col]

if(sum(duplicated(c(var.x.col,var.y.col,df.var.x.col,df.var.y.col,nxi.col,nyi.col))) > 0) warning("Some arguments have common columns in the data set")

if(!is.null(var.x.col) & length(var.x.col) > 1) stop("The var.x.col argument must refer to only one column")
if(!is.null(var.y.col) & length(var.y.col) > 1) stop("The var.y.col argument must refer to only one column")
if(!is.null(nxi.col) & length(nxi.col) > 1) stop("The nxi.col argument must refer to only one column")
if(!is.null(nyi.col) & length(nyi.col) > 1) stop("The nyi.col argument must refer to only one column")
if(!is.null(df.var.x.col) & length(df.var.x.col) > 1) stop("The df.var.x.col argument must refer to only one column")
if(!is.null(df.var.y.col) & length(df.var.y.col) > 1) stop("The df.var.y.col argument must refer to only one column")

if(is.null(var.x.col) & !is.null(df.var.x.col)) warning("df.var.x.col is ignored as var.x.col is not given")
if(is.null(var.y.col) & !is.null(df.var.y.col)) warning("df.var.y.col is ignored as var.y.col is not given")

if(!is.null(var.x.col) & !is.null(var.x.formula)) stop("Choose between var.x.col and var.x.formula (error measurement variance previously calculated in the data set or calculated by a formula)")
if(!is.null(var.y.col) & !is.null(var.y.formula)) stop("Choose between var.y.col and var.y.formula (error measurement variance previously calculated in the data set or calculated by a formula)")
if((is.null(var.x.formula) | is.null(var.y.formula)) & !is.null(xpred)) stop("Not possible to predict outside the X values: use xpred with variances functions")
if(!is.null(var.x.col) & is.null(df.var.x.col)) stop("If the variances are given in a column, the degrees of freedom should also be given (with df.var.x.col)")
if(!is.null(var.y.col) & is.null(df.var.y.col)) stop("If the variances are given in a column, the degrees of freedom should also be given (with df.var.y.col)")
if(!is.null(var.x.formula) & length(grep("X",var.x.formula)) == 0) stop("The variance function (var.x.formula) must contain 'X'")
if(!is.null(var.y.formula) & length(grep("Y",var.y.formula)) == 0) stop("The variance function (var.y.formula) must contain 'Y'")
if((is.null(var.x.col) & is.null(var.x.formula) & length(xcol)==1) | (is.null(var.y.col) & is.null(var.y.formula) & length(ycol)==1)) stop("Not possible to estimate a BLS regression line without any information on the variances (no replicates and no variances in input)")

res=desc.stat(data=data,xcol=xcol,ycol=ycol)
Xij=res$Xij
Yik=res$Yik
x=res$Xi
y=res$Yi
n=res$statistics$N
rr=res$Order.Xi
if(is.null(nxi.col)) nxi=res$nxi else nxi=data[rr,nxi.col]
if(is.null(nyi.col)) nyi=res$nyi else nyi=data[rr,nyi.col]
variances_x=res$variances_x
variances_y=res$variances_y

if(!is.null(nxi.col)) if(!is.numeric(nxi)) stop("Counts of replicates (with nxi.col) must be numeric.")
if(!is.null(nyi.col)) if(!is.numeric(nyi)) stop("Counts of replicates (with nyi.col) must be numeric.")
if(!is.null(df.var.x.col)) if(!is.numeric(data[,df.var.x.col])) stop("Degrees of freedom (with df.var.x.col) must be numeric.")
if(!is.null(df.var.y.col)) if(!is.numeric(data[,df.var.y.col])) stop("Degrees of freedom (with df.var.y.col) must be numeric.")
if(!is.null(var.x.col)) if(!is.numeric(data[,var.x.col])) stop("Measurement error variances (with var.x.col) must be numeric.")
if(!is.null(var.y.col)) if(!is.numeric(data[,var.y.col])) stop("Measurement error variances (with var.y.col) must be numeric.")

if(sum(is.na(nxi) | nxi < 1) > 0) stop("No missing values are allowed in nxi, and all nxi values should be greater than 1.")
if(sum(is.na(nyi) | nyi < 1) > 0) stop("No missing values are allowed in nyi, and all nyi values should be greater than 1.")

if(is.null(var.x.col) & is.null(var.x.formula)) {varX=variances_x/nxi; df.var.x=nxi-1}
if(!is.null(var.x.col) & is.null(var.x.formula)) {variances_x=data[rr,var.x.col]; varX=variances_x/nxi; df.var.x=data[rr,df.var.x.col]}
if(is.null(var.x.col) & !is.null(var.x.formula)) {variances_x=eval(parse(text=sub("X","x",var.x.formula))); varX=variances_x/nxi; df.var.x=rep(n-2,n)}

if(is.null(var.y.col) & is.null(var.y.formula)) {varY=variances_y/nyi; df.var.y=nyi-1}
if(!is.null(var.y.col) & is.null(var.y.formula)) {variances_y=data[rr,var.y.col]; varY=variances_y/nyi; df.var.y=data[rr,df.var.y.col]}
if(is.null(var.y.col) & !is.null(var.y.formula)) {variances_y=eval(parse(text=sub("Y","y",var.y.formula))); varY=variances_y/nyi; df.var.y=rep(n-2,n)}

if(sum(is.na(df.var.x)) > 0) stop("No missing values are allowed in the degrees of freedom")
if(sum(is.na(df.var.y)) > 0) stop("No missing values are allowed in the degrees of freedom")
if(sum(df.var.x <= 0) > 0 | sum(df.var.y <= 0) > 0) warning("Degrees of freedom should be greater than 0")

if(sum(is.na(varX)) > 0) stop("No missing values are allowed in the measurement error variances in X")
if(sum(is.na(varY)) > 0) stop("No missing values are allowed in the measurement error variances in Y")
if(sum(varX < 0) > 0 | sum(varY < 0) > 0)
{
warning("At least one measurement error variance is negative (check the predicted variances). These variances are set to zero.")
varX[which(varX < 0)]=0
varY[which(varY < 0)]=0
}

R=matrix(nrow=2,ncol=2)
g=matrix(nrow=2,ncol=1)
b0=matrix(nrow=2,ncol=1,c(0,1))
diff=matrix(nrow=2,ncol=1)

repeat {
	ei=y-(x*b0[2,1]+b0[1,1])
	ei2=ei^2
	Sei2=varY+b0[2,1]^2*varX
	R1i=1/Sei2
	R2i=x*R1i
	R4i=x^2*R1i	
	R[1,1]=sum(R1i)
	R[1,2]=sum(R2i)
	R[2,1]=sum(R2i)
	R[2,2]=sum(R4i)
	g1i=y*R1i
	g2i=x*y*R1i+b0[2,1]*(ei*R1i)^2*varX
	g[1,1]=sum(g1i)
	g[2,1]=sum(g2i)
	b=solve(R)%*%g
	diff=b0-b
	if ((abs(diff[1,1]) < tol) & (abs(diff[2,1]) < tol))	break()
	b0=b
	#b0=b0-diff/2
	}
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

if(!is.null(var.x.formula) & !is.null(var.y.formula)) xx=c(seq(min(x),max(x),length.out=npoints),xpred) else xx=x
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

if(is.null(var.x.formula) | is.null(var.y.formula))
{
rownames=1:length(x)
} else {
if(is.null(xpred)) rownames=1:npoints else rownames=c(1:npoints,paste("xpred",1:length(xpred)))
}

data.pred=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))
data.pred$X0=xx
data.pred$Ypred=xx.pred

alpha.pred=1-pred.level

data.pred[,CI1]=xx.pred-qt(1-alpha.conf/2,n-2)*sqrt(var.line)
data.pred[,CI2]=xx.pred+qt(1-alpha.conf/2,n-2)*sqrt(var.line)

data.pred[,CB1]=xx.pred-sqrt(qf(conf.level,2,n-2)*Hotelling_correction)*sqrt(var.line)
data.pred[,CB2]=xx.pred+sqrt(qf(conf.level,2,n-2)*Hotelling_correction)*sqrt(var.line)

if(is.null(var.x.col) & !is.null(var.x.formula))
{
var.x.formula=sub("X","xx",var.x.formula)
var_x_single=eval(parse(text=var.x.formula))
df.var.x=rep(n-2,nrow(data.pred))
} else {
var_x_single=varX*nxi
}
if(is.null(var.y.col) & !is.null(var.y.formula))
{
var.y.formula=sub("Y","xx.pred",var.y.formula)
var_y_single=eval(parse(text=var.y.formula))
df.var.y=rep(n-2,nrow(data.pred))
} else {
var_y_single=varY*nyi
}

if(sum(var_x_single < 0) > 0 | sum(var_y_single < 0) > 0)
{
warning("At least one measurement error variance is negative for the predictions. These variances are set to zero.")
var_x_single[which(var_x_single < 0)]=0
var_y_single[which(var_y_single < 0)]=0
}

var_PI=var.line+var_y_single+var_x_single*slope_BLS^2
Numerator=(var.line+var_y_single+slope_BLS^2*var_x_single)^2
Denominator=var.line^2/(n-2)+var_y_single^2/df.var.y+(slope_BLS^2*var_x_single)^2/df.var.x
dfPI=Numerator/Denominator
data.pred[,PI1]=xx.pred-qt(1-alpha.pred/2,dfPI)*sqrt(var_PI)
data.pred[,PI2]=xx.pred+qt(1-alpha.pred/2,dfPI)*sqrt(var_PI)

var_y_q=var_y_single/qy
var_x_q=var_x_single/qx
var_GI=var.line+var_y_q+var_x_q*slope_BLS^2
Numerator=(var.line+var_y_q+slope_BLS^2*var_x_q)^2
Denominator=var.line^2/(n-2)+var_y_q^2/df.var.y+(slope_BLS^2*var_x_q)^2/df.var.x
dfGI=Numerator/Denominator
data.pred[,GI1]=xx.pred-qt(1-alpha.pred/2,dfGI)*sqrt(var_GI)
data.pred[,GI2]=xx.pred+qt(1-alpha.pred/2,dfGI)*sqrt(var_GI)

if(is.null(var.x.formula) | is.null(var.y.formula)) rr=length(x) else rr=npoints

results=list(Xij,Yik,x,y,nxi,nyi,variances_x,variances_y,ell_BLS,Table.BLS,data.pred[1:rr,],data.pred[-(1:rr),])
names(results)=c("Xij","Yik","Xi","Yi","nxi","nyi","variances_x","variances_y","Ellipse.BLS.ht","Estimate.BLS.ht","Pred.BLS.ht","xpred.BLS.ht")
class(results)="BLS.ht"
return(results)

}
