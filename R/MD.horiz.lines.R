#' @export
MD.horiz.lines <-
function(data=NULL,xcol=1,ycol=2,pred.level=0.95,TI.conf.level=0.8)
{

if(length(TI.conf.level) > 1) stop("Only one confidence level is allowed.")
if(length(pred.level) > 1) stop("Only one predictive level is allowed.")
if(TI.conf.level <= 0 | TI.conf.level >= 1) stop("The confidence level must be between 0 and 1 (excluded)")
if(pred.level <= 0 | pred.level >= 1) stop("The predictive level must be between 0 and 1 (excluded)")

res=desc.stat(data=data,xcol=xcol,ycol=ycol)
x=res$Xi
y=res$Yi
n=res$statistics$N
nx=res$statistics$nx
ny=res$statistics$ny
Nx=res$statistics$df.var.x
Ny=res$statistics$df.var.y
varx=res$statistics$'Var X'
vary=res$statistics$'Var Y'
nxprime=1+1/(nx*n)
nyprime=1+1/(ny*n)

if((nx == 1 & ny > 1) | (nx > 1 & ny == 1)) stop("Measurements of both methods should be unreplicated or replicated")

data=cbind(x,y)
m=apply(data,1,mean)
d=data[,2]-data[,1]
data.MD=as.data.frame(cbind(m,d))
colnames(data.MD)=c("Means","Differences")

Dmean=mean(d)
varD=var(d)
if(nx == 1 & ny == 1) varD.single=varD else varD.single=varx+vary

rownames=paste("Differences")
colnames=c("Mean","Std Dev","Std Dev single","Min","Q1","Q2","Q3","Max","N")
Table.Differences=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))
Table.Differences$Mean=Dmean
Table.Differences$'Std Dev'=sqrt(varD)
Table.Differences$'Std Dev single'=sqrt(varD.single)
Table.Differences$'Min'=min(d,na.rm=T)
Table.Differences$'Q1'=quantile(d,0.25,na.rm=T)
Table.Differences$'Q2'=quantile(d,0.5,na.rm=T)
Table.Differences$'Q3'=quantile(d,0.75,na.rm=T)
Table.Differences$'Max'=max(d,na.rm=T)
Table.Differences$'N'=nrow(data)

rownames=c("AI","b TI","bg TI")
colnames=c("Mean differences","Lower","Upper","Meaning","Remark")
Intervals=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))

alpha.pred=1-pred.level

if(nx==1 & ny==1)
{
	AI1=Dmean-qnorm(1-alpha.pred/2)*sqrt(varD)
	AI2=Dmean+qnorm(1-alpha.pred/2)*sqrt(varD)

	bTI1=Dmean-qt(1-alpha.pred/2,n-1)*sqrt(varD)*sqrt(1+1/n)
	bTI2=Dmean+qt(1-alpha.pred/2,n-1)*sqrt(varD)*sqrt(1+1/n)

	bgTI1=Dmean-qnorm(1-alpha.pred/2)*sqrt(varD)*sqrt(1+1/n)*sqrt((n-1)/qchisq(1-TI.conf.level,n-1))
	bgTI2=Dmean+qnorm(1-alpha.pred/2)*sqrt(varD)*sqrt(1+1/n)*sqrt((n-1)/qchisq(1-TI.conf.level,n-1))

} else {

	varDs=varD+(1-1/nx)*varx+(1-1/ny)*vary
	AI1=Dmean-qnorm(1-alpha.pred/2)*sqrt(varDs)
	AI2=Dmean+qnorm(1-alpha.pred/2)*sqrt(varDs)

	r.omega=df.WS(variances=c(varx,vary),k=c(nxprime,nyprime),dfs=c(Nx,Ny))
	bTI1=Dmean-qt(1-alpha.pred/2,r.omega)*sqrt(vary*nyprime+varx*nxprime)
	bTI2=Dmean+qt(1-alpha.pred/2,r.omega)*sqrt(vary*nyprime+varx*nxprime)

	H1=Nx/qchisq(1-TI.conf.level,Nx)-1
	H2=Ny/qchisq(1-TI.conf.level,Ny)-1
	H_all=nxprime*H1^2*varx^2+nyprime*H2^2*vary^2
	bgTI1=Dmean-qnorm(1-alpha.pred/2)*sqrt(vary*nyprime+varx*nxprime)*sqrt(1+(1/varD.single)*sqrt(H_all))
	bgTI2=Dmean+qnorm(1-alpha.pred/2)*sqrt(vary*nyprime+varx*nxprime)*sqrt(1+(1/varD.single)*sqrt(H_all))
}

Intervals$'Mean differences'=Dmean
Intervals$Lower=c(AI1,bTI1,bgTI1)
Intervals$Upper=c(AI2,bTI2,bgTI2)
Intervals$Meaning=c(paste(100*pred.level,"% Agreement interval",sep=""),paste(100*pred.level,"% Tolerance interval",sep=""),paste(100*pred.level,"% Tolerance interval (with ",100*TI.conf.level,"% confidence level)",sep=""))
texte1="Too narrow, especially for small sample sizes. Not robust to points outliyng the distribution of the differences (with replicated data)"
texte2="Exact solution (unreplicated data). Robust (with replicates)"
texte3="Add a confidence level to the beta expectation TI. Robust (with replicates)"
Intervals$Remark=c(texte1,texte2,texte3)

results=list(data.MD,Table.Differences,Intervals)
names(results)=c("data.MD","Table.Differences","Intervals.horiz.lines")
class(results)="MD.horiz.lines"
return(results)
}
