#' @export
FullCIs.MD <-
function(data=NULL,xcol=1,ycol=2,conf.level=0.95,npoints=1000,nlambdas=13)
{

if(length(conf.level) > 1) stop("Only one confidence level is allowed.")
if(conf.level <= 0 | conf.level >= 1) stop("The confidence level must be between 0 and 1 (excluded)")
if(length(npoints) > 1) stop("npoints must be an integer (at least 10)")
if(floor(npoints) != ceiling(npoints) | npoints < 10 | !is.numeric(npoints)) stop("npoints must be an integer (at least 10)")
if(length(nlambdas) > 1) stop("nlambdas must be an integer")
if(floor(nlambdas) != ceiling(nlambdas) | nlambdas < 2 | !is.numeric(nlambdas)) stop("nlambdas must be an integer (at least 2)")

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

d=y-x
m=(x+y)/2
y=d
x=m
lambdaxy=4
xmean=mean(x)
xcentre=x-xmean
Sxx=sum((x-xmean)^2)
ymean=mean(y)
ycentre=y-ymean
Syy=sum((y-ymean)^2)
Sxy=sum(xcentre*ycentre)

alpha=1-conf.level
Hotelling_correction=(2*(n-1))/(n-2)

if(nlambdas == 13) epsilons=c(0.000001,0.1,0.175,0.25,0.375,0.5,0.625,0.75,0.875,1,1.25,1.5,1.75,2,200000) else epsilons=c(0.000001,10^seq(log10(0.1),log10(2),length.out=nlambdas),200000)
deltas=rev(epsilons)
lambdas=epsilons/deltas
CBLS.corrs=(lambdas-1)/(lambdas+1)

ell_CBLS_CB_all=array(dim=c(npoints,2,length(epsilons)),dimnames=list(1:npoints,c("Slopes","Intercepts"),paste("Lambda",1:length(epsilons),sep="")))
slopes_all=matrix(nrow=length(epsilons),ncol=4) # 7: Estimate, CI lower, CI upper, PVal
intercepts_all=matrix(nrow=length(epsilons),ncol=4) # 4: Estimate, CI lower, CI upper, PVal
joints_all=matrix(nrow=length(epsilons),ncol=1) # 1: PVal

xx=c(seq(min(x),max(x),length.out=npoints)) # To smooth the hyperbolic intervals
CI1=paste(conf.level*100,"% CI Lower",sep="")
CI2=paste(conf.level*100,"% CI Upper",sep="")
CB1=paste(conf.level*100,"% CB Lower",sep="")
CB2=paste(conf.level*100,"% CB Upper",sep="")

colnames=c("X0","Ypred",CI1,CI2,CB1,CB2)
hyperbola_CBLS_all=array(dim=c(npoints,length(colnames),length(epsilons)),dimnames=list(1:npoints,colnames,paste("Lambda",1:length(epsilons),sep="")))
hyperbola_CBLS_all[,"X0",]=xx

for (j in 1:length(deltas))
{
	CBLS.corr=CBLS.corrs[j]
	
	slopes_all[j,1]=(Syy-4*Sxx+sqrt((Syy-4*Sxx)^2-4*(Sxy-2*CBLS.corr*Sxx)*(2*CBLS.corr*Syy-4*Sxy)))/(2*(Sxy-2*CBLS.corr*Sxx))
	intercepts_all[j,1]=ymean-xmean*slopes_all[j,1]

	varX=rep(1,n)
	varY=rep(4,n)
	CBLS.cov=CBLS.corr*sqrt(varX[1])*sqrt(varY[1])

	R=matrix(nrow=2,ncol=2)
	g=matrix(nrow=2,ncol=1)
	#diff=matrix(nrow=2,ncol=1)
	#repeat {
		ei=y-(x*slopes_all[j,1]+intercepts_all[j,1])
	#	ei2=ei^2
		Sei2=varY+slopes_all[j,1]^2*varX-2*slopes_all[j,1]*CBLS.cov
		R1i=1/Sei2
		R2i=x*R1i
		R4i=x^2*R1i	
		R[1,1]=sum(R1i)
		R[1,2]=sum(R2i)
		R[2,1]=sum(R2i)
		R[2,2]=sum(R4i)
		g1i=y*R1i
		g2i=x*y*R1i+(ei*R1i)^2*(slopes_all[j,1]*varX-CBLS.cov)
		g[1,1]=sum(g1i)
		g[2,1]=sum(g2i)
		b=solve(R)%*%g
	#	diff=b0-b
	#	b0=b0-diff/2
	#	if ((abs(diff[1,1])<0.000001) & (abs(diff[2,1])<0.000001))
	#	break()
	#	}
	#slopes_all[j,4]=b[2,1]
	#intercepts_all[j,1]=b[1,1]
	s2i=sum(((y-(x*slopes_all[j,1]+intercepts_all[j,1]))^2)/Sei2)/(n-2)
	varcovari=solve(R)*s2i
	S_slope_CBLS_diagram=sqrt(varcovari[2,2])
	S_intercept_CBLS_diagram=sqrt(varcovari[1,1])
	covar_slope_intercept_CBLS_diagram=varcovari[1,2]
	cov_matrix_CBLS_diagram=matrix(nrow=2,ncol=2,c(S_slope_CBLS_diagram^2,covar_slope_intercept_CBLS_diagram,covar_slope_intercept_CBLS_diagram,S_intercept_CBLS_diagram^2))
	slopes_all[j,2]=slopes_all[j,1]-qt(1-alpha/2,n-2)*S_slope_CBLS_diagram
	slopes_all[j,3]=slopes_all[j,1]+qt(1-alpha/2,n-2)*S_slope_CBLS_diagram
	slopes_all[j,4]=(1-pt(abs(slopes_all[j,1]-0)/S_slope_CBLS_diagram,n-2))*2
	intercepts_all[j,2]=intercepts_all[j,1]-qt(1-alpha/2,n-2)*S_intercept_CBLS_diagram
	intercepts_all[j,3]=intercepts_all[j,1]+qt(1-alpha/2,n-2)*S_intercept_CBLS_diagram
	intercepts_all[j,4]=(1-pt(abs(intercepts_all[j,1])/S_intercept_CBLS_diagram,n-2))*2
	ell_CBLS_CB_all[,,j]=ellipse(cov_matrix_CBLS_diagram,centre=c(slopes_all[j,1],intercepts_all[j,1]),npoints=npoints, t = sqrt(qf(conf.level,2,n-2)*Hotelling_correction))

	F1=matrix(nrow=1,ncol=2,c(slopes_all[j,1]-0,intercepts_all[j,1]))
	F2=solve(cov_matrix_CBLS_diagram)
	F3=t(F1)
	Fell=F1%*%F2%*%F3
	pval_ell_CBLS=1-pf(F1%*%F2%*%F3/Hotelling_correction,2,n-2)
	joints_all[j,1]=pval_ell_CBLS

	hyperbola_CBLS_all[,"Ypred",j]=xx*slopes_all[j,1]+intercepts_all[j,1]
	var.line=S_intercept_CBLS_diagram^2+S_slope_CBLS_diagram^2*xx^2+2*xx*covar_slope_intercept_CBLS_diagram

	hyperbola_CBLS_all[,CI1,j]=hyperbola_CBLS_all[,"Ypred",j]-qt(1-alpha/2,n-2)*sqrt(var.line)
	hyperbola_CBLS_all[,CI2,j]=hyperbola_CBLS_all[,"Ypred",j]+qt(1-alpha/2,n-2)*sqrt(var.line)

	hyperbola_CBLS_all[,CB1,j]=hyperbola_CBLS_all[,"Ypred",j]-sqrt(qf(conf.level,2,n-2)*Hotelling_correction)*sqrt(var.line)
	hyperbola_CBLS_all[,CB2,j]=hyperbola_CBLS_all[,"Ypred",j]+sqrt(qf(conf.level,2,n-2)*Hotelling_correction)*sqrt(var.line)

}

lambdas[1]=0
lambdas[length(lambdas)]=Inf

slopes=cbind(lambdas,(lambdas-1)/(lambdas+1),slopes_all)
colnames(slopes)=c("Ratio variances XY","Correlation MD","Estimate",paste("Lower ",conf.level*100,"%CI",sep=""),paste("Upper ",conf.level*100,"%CI",sep=""),"pvalue")
intercepts=cbind(lambdas,(lambdas-1)/(lambdas+1),intercepts_all)
colnames(intercepts)=c("Ratio variances XY","Correlation MD","Estimate",paste("Lower ",conf.level*100,"%CI",sep=""),paste("Upper ",conf.level*100,"%CI",sep=""),"pvalue")
joints=cbind(lambdas,(lambdas-1)/(lambdas+1),joints_all)
colnames(joints)=c("Ratio variances XY","Correlation MD","pvalue")

rr=which(slopes[,"Ratio variances XY"] == Inf)
slopes[rr,"Correlation MD"]=1
intercepts[rr,"Correlation MD"]=1
joints[rr,"Correlation MD"]=1

res=desc.stat(data=data,xcol=xcol,ycol=ycol)
x=res$Xi
y=res$Yi
d=y-x
m=(x+y)/2

data.means=as.data.frame(cbind(x,y,m,d))
colnames(data.means)=c("X","Y","M","D")

results=list(data.means,ell_CBLS_CB_all,as.data.frame(slopes),as.data.frame(intercepts),as.data.frame(joints),hyperbola_CBLS_all)
names(results)=c("Data.means","Ellipses.CB","Slopes","Intercepts","Joints","Hyperbolic.intervals")
class(results)="CIs.MD"
return(results)
}
