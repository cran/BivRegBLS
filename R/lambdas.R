#' @export
lambdas <-
function(data=NULL,xcol=NULL,ycol=NULL,conf.level=0.95)
{

if (length(xcol) <= 1 | length(ycol) <= 1 | is.null(xcol) | is.null(ycol)) stop("No replicates found. xcol and ycol must be, each, at least of length two")
if(sum(conf.level <= 0 | conf.level >= 1) > 0) stop("The confidence levels must be between 0 and 1 (excluded)")
if(length(conf.level) > 1) warning("At least 2 confidence levels are specified. Choose only one level if these results need to be plugged into GraphFullCIs.XY or GraphFullCIs.MD functions")

res=desc.stat(data=data,xcol=xcol,ycol=ycol)
Nx=res$statistics$df.var.x
Ny=res$statistics$df.var.y
nx=res$statistics$nx
ny=res$statistics$ny
vary=res$statistics$'Var Y'
varx=res$statistics$'Var X'
k_unb=(Nx-2)/Nx

rownames=c("Lambda","LambdaXY")
name_CI1=paste("Lower ",conf.level*100,"%CI",sep="")
name_CI2=paste("Upper ",conf.level*100,"%CI",sep="")
colnames=c("Estimate","Unbiased estimate",name_CI1,name_CI2,"pvalue")
lambdas.Table=as.data.frame(matrix(nrow=length(rownames),ncol=length(colnames),dimnames=list(rownames,colnames)))

lambdas.Table$Estimate[1]=vary/varx
lambdas.Table$'Unbiased estimate'[1]=lambdas.Table$Estimate[1]*k_unb

alpha=1-conf.level
lambdas.Table[1,name_CI1]=lambdas.Table$Estimate[1]/qf(1-alpha/2,Ny,Nx)
lambdas.Table[1,name_CI2]=lambdas.Table$Estimate[1]/qf(alpha/2,Ny,Nx)
lambdas.Table$pvalue[1]=min(pf(lambdas.Table$Estimate[1],Ny,Nx)*2,(1-pf(lambdas.Table$Estimate[1],Ny,Nx))*2)

k=nx/ny
lambdas.Table$Estimate[2]=lambdas.Table$Estimate[1]*k
lambdas.Table[2,name_CI1]=lambdas.Table[1,name_CI1]*k
lambdas.Table[2,name_CI2]=lambdas.Table[1,name_CI2]*k
lambdas.Table$'Unbiased estimate'[2]=lambdas.Table$Estimate[2]*k_unb
lambdas.Table$pvalue[2]=min(pf(lambdas.Table$Estimate[1]*k,Ny,Nx)*2,(1-pf(lambdas.Table$Estimate[1]*k,Ny,Nx))*2)

results=list(lambdas.Table)
names(results)=c("lambdas.Table")
class(results)="lambdas"
return(results)
}
