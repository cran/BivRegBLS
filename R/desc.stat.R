#' @importFrom stats cor pf pt qchisq qf qnorm qt quantile var
#' @export
desc.stat <-
function(data=NULL,xcol=1,ycol=2,IDcol=NULL)
{
if(is.null(data)) stop("No data set is specified by the argument 'data'")
if(is.null(colnames(data))) colnames(data)=paste0("V",1:ncol(data))
if(sum((is.element(xcol,colnames(data)) | xcol <= ncol(data)) == F) > 0) stop ("The X columns (xcol) are not in the data set")
if(sum((is.element(ycol,colnames(data)) | ycol <= ncol(data)) == F) > 0) stop ("The Y columns (ycol) are not in the data set")

if(is.numeric(xcol)) xcol=colnames(data)[xcol]
if(is.numeric(ycol)) ycol=colnames(data)[ycol]
if(!is.null(IDcol) & is.numeric(IDcol)) IDcol=colnames(data)[IDcol]

if(is.null(IDcol))
{
if(sum(xcol %in% ycol) > 0) stop("X and Y columns must be separate (no common columns)")
} else {
if(sum((is.element(IDcol,colnames(data)) | IDcol <= ncol(data)) == F) > 0) stop ("The ID column (IDcol) is not in the data set")
if(sum(xcol %in% ycol, xcol %in% IDcol, ycol %in% IDcol) > 0) stop("X, Y and ID columns must be separate (no common columns)")
if(length(IDcol) > 1) stop("The ID column must refer to only one column")
}

xi=as.matrix(data[,xcol])
yi=as.matrix(data[,ycol])

if(!is.numeric(xi) | !is.numeric(yi)) stop("Data are not numeric")

x=apply(xi,1,mean,na.rm=TRUE)
nxi=apply(!is.na(xi),1,sum)
variances_x=apply(xi,1,var,na.rm=TRUE)

y=apply(yi,1,mean,na.rm=TRUE)
nyi=apply(!is.na(yi),1,sum)
variances_y=apply(yi,1,var,na.rm=TRUE)

data.summary=as.data.frame(cbind(x,y))

if(is.null(IDcol)) data.summary$IDs=rep(1,nrow(data.summary)) else data.summary$IDs=data[,IDcol]

data.summary$nxi=nxi
data.summary$nyi=nyi
data.summary$variances_x=variances_x
data.summary$variances_y=variances_y

ww=which(is.na(data.summary$x) | is.na(data.summary$y))
if(length(ww) > 0)
{
warning(paste(length(ww)," row(s) removed because of missingness in mean measures",sep=""))
xi=as.matrix(xi[-ww,])
yi=as.matrix(yi[-ww,])
data.summary=data.summary[-ww,]
}

if (min(data.summary$nxi) != max(data.summary$nxi)) warning("Count of replicates X not constant, the maximum is considered")
if (min(data.summary$nyi) != max(data.summary$nyi)) warning("Count of replicates Y not constant, the maximum is considered")

if(!is.null(IDcol)) Nb.IDs=length(levels(as.factor(data.summary$IDs))) else Nb.IDs=1
if(Nb.IDs > 30) stop("Not more than 30 IDs are allowed")

rr=order(data.summary$x)
data.summary=data.summary[rr,]
xi=xi[rr,]
yi=yi[rr,]

colnames=c("N","nx","ny","df.var.x","df.var.y","Xmean","Sxx","Var X","Ymean","Syy","Sxy","Var Y","R","R2","Min.X","Q1.X","Q2.X","Q3.X","Max.X","Min.Y","Q1.Y","Q2.Y","Q3.Y","Max.Y")

if(Nb.IDs == 1)
{
stats=as.data.frame(matrix(nrow=1,ncol=length(colnames),dimnames=list("Statistics",colnames)))
IDs=rep(1,nrow(data))
} else {
stats=as.data.frame(matrix(nrow=Nb.IDs,ncol=length(colnames),dimnames=list(paste("Group",levels(as.factor(data.summary$IDs))),colnames)))
}

for (i in 1:Nb.IDs)
{
subdata=which(data.summary$IDs==levels(as.factor(data.summary$IDs))[i])
nx=max(data.summary$nxi[subdata],na.rm=T)
var_x=sum((data.summary$nxi[subdata]-1)*data.summary$variances_x[subdata],na.rm=T)/sum(data.summary$nxi[subdata]-1)
ny=max(data.summary$nyi[subdata],na.rm=T)
var_y=sum((data.summary$nyi[subdata]-1)*data.summary$variances_y[subdata],na.rm=T)/sum(data.summary$nyi[subdata]-1)

N=length(subdata)
df.var.x=sum(data.summary$nxi[subdata]-1)
df.var.y=sum(data.summary$nyi[subdata]-1)

stats$N[i]=N
stats$nx[i]=nx
stats$ny[i]=ny
stats$df.var.x[i]=df.var.x
stats$df.var.y[i]=df.var.y
stats$Xmean[i]=mean(data.summary$x[subdata])
stats$'Var X'[i]=var_x
stats$Sxx[i]=sum((data.summary$x[subdata]-stats$Xmean[i])^2)
stats$Ymean[i]=mean(data.summary$y[subdata])
stats$Syy[i]=sum((data.summary$y[subdata]-stats$Ymean[i])^2)
stats$'Var Y'[i]=var_y
stats$Sxy[i]=sum((data.summary$x[subdata]-stats$Xmean[i])*(data.summary$y[subdata]-stats$Ymean[i]))
stats$R[i]=cor(data.summary$x[subdata],data.summary$y[subdata])
stats$R2[i]=stats$R[i]^2
stats$Min.X[i]=min(data.summary$x[subdata])
stats$Q1.X[i]=quantile(data.summary$x[subdata],0.25)
stats$Q2.X[i]=quantile(data.summary$x[subdata],0.5)
stats$Q3.X[i]=quantile(data.summary$x[subdata],0.75)
stats$Max.X[i]=max(data.summary$x[subdata])
stats$Min.Y[i]=min(data.summary$y[subdata])
stats$Q1.Y[i]=quantile(data.summary$y[subdata],0.25)
stats$Q2.Y[i]=quantile(data.summary$y[subdata],0.5)
stats$Q3.Y[i]=quantile(data.summary$y[subdata],0.75)
stats$Max.Y[i]=max(data.summary$y[subdata])
}# for i

results=list(xi,yi,data.summary$x,data.summary$y,data.summary$IDs,data.summary$nxi,data.summary$nyi,data.summary$variances_x,data.summary$variances_y,rr,stats)
names(results)=c("Xij","Yik","Xi","Yi","IDs","nxi","nyi","variances_x","variances_y","Order.Xi","statistics")
class(results)="desc.stat"
return(results)
}
