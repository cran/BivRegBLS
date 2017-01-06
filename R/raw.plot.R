#' @export
raw.plot <-
function(data.plot=NULL,xname="X",yname="Y",graph="XY.means",col.ID=NULL,pch.ID=NULL,...)
{

if(class(data.plot) != "desc.stat") stop("data.plot must be a desc.stat class object (obtained after the desc.stat function)")
if(length(xname) > 1 | !is.character(xname)) stop("xname must be a character name")
if(length(yname) > 1 | !is.character(yname)) stop("yname must be a character name")
if(length(graph) > 1) stop("Only one type of graph is allowed in 'graph'")
if(!(graph %in% c("XY.means", "XY.points", "XY.bar.range", "XY.bar.SEM", "XY.bar.SD","MD.means"))) stop("Only the 6 following options are allowed in graph: XY.means, XY.points, XY.bar.range, XY.bar.SEM, XY.bar.SD")

Xij=as.matrix(data.plot$Xij)
Yik=as.matrix(data.plot$Yik)
Xi=data.plot$Xi
Yi=data.plot$Yi
IDs=data.plot$IDs
nxi=data.plot$nxi
nyi=data.plot$nyi
variances_x=data.plot$variances_x
variances_y=data.plot$variances_y
N=nrow(Xij)
nx=max(nxi,na.rm=T)
ny=max(nyi,na.rm=T)
Xi.min=apply(Xij,1,min,na.rm=T)
Xi.max=apply(Xij,1,max,na.rm=T)
Yi.min=apply(Yik,1,min,na.rm=T)
Yi.max=apply(Yik,1,max,na.rm=T)
Xi.sd1=Xi-sqrt(variances_x)
Xi.sd2=Xi+sqrt(variances_x)
Yi.sd1=Yi-sqrt(variances_y)
Yi.sd2=Yi+sqrt(variances_y)
Xi.sem1=Xi-sqrt(variances_x/nxi)
Xi.sem2=Xi+sqrt(variances_x/nxi)
Yi.sem1=Yi-sqrt(variances_y/nyi)
Yi.sem2=Yi+sqrt(variances_y/nyi)
Mi=(Xi+Yi)/2
Di=Yi-Xi

Nb.IDs=length(levels(as.factor(IDs)))
if(is.null(pch.ID)) pch.ID=1:Nb.IDs else pch.ID=rep(pch.ID,Nb.IDs)
if(is.null(col.ID)) col.ID=rep(1,Nb.IDs) else col.ID=rep(col.ID,Nb.IDs)

if(sum(pch.ID > 25) >= 1) stop("Don't use pch value greater than 25. Change the pch.ID value")
if(nx == 1 & ny == 1 & graph %in% c("XY.bar.range", "XY.bar.SEM", "XY.bar.SD")) stop("This kind of graph is not suitable for unreplicated data")

arguments <- list(...)
names.arguments=names(arguments)

graph.raw=function()
{
	if(graph != "MD.means")
	{
		if(!("xlim" %in% names.arguments))
		{
		min_x=min(Xij,Yik,na.rm=T)
		max_x=max(Xij,Yik,na.rm=T)
		} else {
		min_x=arguments$xlim[1]
		max_x=arguments$xlim[2]
		}
		if(!("ylim" %in% names.arguments))
		{
		min_y=min(Xij,Yik,na.rm=T)
		max_y=max(Xij,Yik,na.rm=T)
		} else {
		min_y=arguments$ylim[1]
		max_y=arguments$ylim[2]
		}
		if(!("xlab" %in% names.arguments)) xlab="" else xlab=arguments$xlab
		if(!("ylab" %in% names.arguments)) ylab="" else ylab=arguments$ylab
		arguments.new=arguments
		arguments.new[c("xlim","ylim","xlab","ylab")]=NULL

	if(Nb.IDs > 1) par(mar=c(6,6,6,6))
	do.call(plot,args=c(list(x=1:10),list(y=1:10),list(type="n"),list(xlim=c(min_x,max_x)),list(ylim=c(min_y,max_y)),list(xlab=xlab),list(ylab=ylab),arguments.new))
	abline(a=0,b=1,lty=2)

	if(graph == "XY.means")
	{
		if(!("xlab" %in% names.arguments)) title(xlab=xname)
		if(!("ylab" %in% names.arguments)) title(ylab=yname)
		for(i in 1:Nb.IDs)
		{
		subdata=which(IDs==levels(as.factor(IDs))[i])
		points(Xi[subdata],Yi[subdata],pch=pch.ID[i],col=col.ID[i])
		}
	}
	if(graph == "XY.points")
	{
		if(!("xlab" %in% names.arguments)) title(xlab=xname)
		if(!("ylab" %in% names.arguments)) title(ylab=yname)
		for(i in 1:Nb.IDs)
		{
		subdata=which(IDs==levels(as.factor(IDs))[i])
		for(j in 1:nx) points(Xij[subdata,j],Yi[subdata],pch=pch.ID[i],col=col.ID[i])
		for(k in 1:ny) points(Xi[subdata],Yik[subdata,k],pch=pch.ID[i],col=col.ID[i])
		}
	}
	if(graph == "XY.bar.range")
	{
		if(!("xlab" %in% names.arguments)) title(xlab=paste(xname," (Min-Max)",sep=""))
		if(!("ylab" %in% names.arguments)) title(ylab=paste(yname," (Min-Max)",sep=""))
		for(i in 1:Nb.IDs)
		{
		subdata=which(IDs==levels(as.factor(IDs))[i])
		for(s in 1:length(subdata)) segments(Xi.min[subdata][s],Yi[subdata][s],Xi.max[subdata][s],Yi[subdata][s],col=col.ID[i])
		for(s in 1:length(subdata)) segments(Xi[subdata][s],Yi.min[subdata][s],Xi[subdata][s],Yi.max[subdata][s],col=col.ID[i])
		}
	}
	if(graph == "XY.bar.SD")
	{
		if(!("xlab" %in% names.arguments)) title(xlab=paste(xname," (+/- SD)",sep=""))
		if(!("ylab" %in% names.arguments)) title(ylab=paste(yname," (+/- SD)",sep=""))
		for(i in 1:Nb.IDs)
		{
		subdata=which(IDs==levels(as.factor(IDs))[i])
		for(s in 1:length(subdata)) segments(Xi.sd1[subdata][s],Yi[subdata][s],Xi.sd2[subdata][s],Yi[subdata][s],col=col.ID[i])
		for(s in 1:length(subdata)) segments(Xi[subdata][s],Yi.sd1[subdata][s],Xi[subdata][s],Yi.sd2[subdata][s],col=col.ID[i])
		}
	}
	if(graph == "XY.bar.SEM")
	{
		if(!("xlab" %in% names.arguments)) title(xlab=paste(xname," (+/- SEM)",sep=""))
		if(!("ylab" %in% names.arguments)) title(ylab=paste(yname," (+/- SEM)",sep=""))
		for(i in 1:Nb.IDs)
		{
		subdata=which(IDs==levels(as.factor(IDs))[i])
		for(s in 1:length(subdata)) segments(Xi.sem1[subdata][s],Yi[subdata][s],Xi.sem2[subdata][s],Yi[subdata][s],col=col.ID[i])
		for(s in 1:length(subdata)) segments(Xi[subdata][s],Yi.sem1[subdata][s],Xi[subdata][s],Yi.sem2[subdata][s],col=col.ID[i])
		}
	}
	legend("bottomright",lty=2,c("Y=X"),bty="n")
	if(Nb.IDs > 1) legend(max_x+(max_x-min_x)*0.05,max_y+(max_y-min_y)*0.13,levels(as.factor(IDs)),xpd=NA,pch=pch.ID,col=col.ID,horiz=F,bty="n")
	}#if
	
	if(graph == "MD.means")
	{
	if(!("ylab" %in% names.arguments)) y.title=paste("Differences (D = ",yname," - ",xname,")",sep="") else y.title=arguments$ylab
	if(!("xlab" %in% names.arguments)) x.title=paste("Averages (M = (",xname," + ",yname,")/2)",sep="") else x.title=arguments$xlab
	if(!("xlim" %in% names.arguments))
	{
	min_x=min(Mi,na.rm=T)
	max_x=max(Mi,na.rm=T)
	} else {
	min_x=arguments$xlim[1]
	max_x=arguments$xlim[2]
	}
	if(!("ylim" %in% names.arguments))
	{
	min_y=min(Di,na.rm=T)
	max_y=max(Di,na.rm=T)
	leg=(max_y-min_y)*0.03
	min_y=min_y-leg
	} else {
	min_y=arguments$ylim[1]
	max_y=arguments$ylim[2]
	}

	arguments.new=arguments
	arguments.new[c("xlim","ylim","xlab","ylab")]=NULL

	do.call(plot,args=c(list(x=1:10),list(y=1:10),list(type="n"),list(xlim=c(min_x,max_x)),list(ylim=c(min_y,max_y)),list(xlab=x.title),list(ylab=y.title),arguments.new))
	abline(a=0,b=0,lty=2)
	
		for(i in 1:Nb.IDs)
		{
		subdata=which(IDs==levels(as.factor(IDs))[i])
		points(Mi[subdata],Di[subdata],pch=pch.ID[i],col=col.ID[i])
		}
	legend("bottomright",lty=2,c("D=0"),bty="n")
	}
	
}#graph.raw

graph.raw()

}
