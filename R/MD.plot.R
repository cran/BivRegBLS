#' @export
MD.plot <-
function(results=NULL,xname="X",yname="Y",antilog=NULL,accept.int=0,accept.int.perc=FALSE,accept.int.threshold=NULL,include.int=TRUE,graph.int=c("CB","PI"),graph.horiz.int=c("bTI","bgTI"),col.CBLS=1,col.CI=2,col.CB=3,col.PI=4,col.GI=5,col.bTI=3,col.bgTI=4,lty.CBLS=1,lty.CI=1,lty.CB=1,lty.PI=1,lty.GI=1,lty.bTI=1,lty.bgTI=1,...)
{

if(!is.null(antilog)) {if(antilog != 10 & antilog != "e") stop('antilog must be equal to "e" or 10')}
if(class(results) != "CBLS" & class(results) != "MD.horiz.lines") stop("results must be a CBLS or MD.horiz.lines class object (obtained after the CBLS or MD.horiz.lines functions)")
if(length(xname) > 1 | !is.character(xname)) stop("xname must be a character name")
if(length(yname) > 1 | !is.character(yname)) stop("yname must be a character name")
if( sum(graph.int != "CI" & graph.int != "CB" & graph.int != "PI" & graph.int != "GI") > 0) stop("Only CI, CB, PI or GI are allowed in graph.int")

if(length(accept.int) > 2) stop("maximum two values are allowed in accept.int")
if(!is.numeric(accept.int)) stop("accept.int must be numeric")
if(length(accept.int.perc) > 1) stop("accept.int.perc must be TRUE or FALSE")
if(!(accept.int.perc %in% c(TRUE,FALSE))) stop("accept.int.perc must be TRUE or FALSE")
if(length(include.int) > 1) stop("include.int must be TRUE or FALSE")
if(!(include.int %in% c(TRUE,FALSE))) stop("include.int must be TRUE or FALSE")
if(!is.null(accept.int.threshold)) {if(length(accept.int.threshold) > 1 | !is.numeric(accept.int.threshold)) stop("accept.int.threshold must be one numeric value")}
if(length(accept.int) == 2 & is.null(accept.int.threshold)) stop("If 2 values are used in the acceptance interval, then a threshold in X must be given by accept.int.threshold")

if(!is.null(antilog)) warning("accept.int.perc is ignored if an antilog is applied")

if(is.null(antilog))
{
y.title=paste("Differences: D = ",yname," - ",xname,sep="")
x.title=paste("Averages: M = (",xname," + ",yname,")/2",sep="")
leg0="D=0"
} else {
y.title=paste("Ratios: R = ",yname," / ",xname,sep="")
x.title=paste("Geom. Means: {",xname,", ",yname,"}",sep="")
leg0="R=1"
}

arguments <- list(...)
names.arguments=names(arguments)

if(class(results) == "CBLS")
{
res=results
if(!is.null(antilog)) res=antilog.pred(results=res,base=antilog)
x=res$Mi
y=res$Di
Pred.CBLS=res$Pred.CBLS
if(!is.null(accept.int.threshold)) {if(accept.int.threshold < min(x) | accept.int.threshold > max(x)) stop("accept.int.threshold should be between the minimum and maximum of (X+Y)/2 values, or use only one acceptance interval value.")}

coltitle=colnames(Pred.CBLS)
conf.level=substring(coltitle[grep("CI",coltitle)[1]],1,regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-1)
pred.level=substring(coltitle[grep("PI",coltitle)[1]],1,regexpr("%",coltitle[grep("PI",coltitle)[1]])[[1]]-1)
CI1=paste(conf.level,"% CI Lower",sep="")
CI2=paste(conf.level,"% CI Upper",sep="")
CB1=paste(conf.level,"% CB Lower",sep="")
CB2=paste(conf.level,"% CB Upper",sep="")
PI1=paste(pred.level,"% PI Lower",sep="")
PI2=paste(pred.level,"% PI Upper",sep="")
GI1=paste(pred.level,"% GI Lower",sep="")
GI2=paste(pred.level,"% GI Upper",sep="")

graph.MD=function()
{

if(!is.null(antilog)) diff.H0=1 else diff.H0=0

if(length(accept.int) == 1 & !is.null(antilog))
{
x.polygon=c(min(x),min(x),max(x),max(x))
y.polygon=c(diff.H0-accept.int/100,1/(diff.H0-accept.int/100),1/(diff.H0-accept.int/100),diff.H0-accept.int/100)
}
if(accept.int.perc == FALSE & length(accept.int) == 1 & is.null(antilog))
{
x.polygon=c(min(x),min(x),max(x),max(x))
y.polygon=c(diff.H0-accept.int,diff.H0+accept.int,diff.H0+accept.int,diff.H0-accept.int)
}
if(length(accept.int) == 2 & !is.null(antilog))
{
x.polygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
y.polygon=c(diff.H0-accept.int[1]/100,1/(diff.H0-accept.int[1]/100),1/(diff.H0-accept.int[1]/100),1/(diff.H0-accept.int[2]/100),1/(diff.H0-accept.int[2]/100),diff.H0-accept.int[2]/100,diff.H0-accept.int[2]/100,diff.H0-accept.int[1]/100)
}
if(accept.int.perc == FALSE & length(accept.int) == 2 & is.null(antilog))
{
x.polygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
y.polygon=c(diff.H0-accept.int[1],diff.H0+accept.int[1],diff.H0+accept.int[1],diff.H0+accept.int[2],diff.H0+accept.int[2],diff.H0-accept.int[2],diff.H0-accept.int[2],diff.H0-accept.int[1])
}
if(accept.int.perc == TRUE & length(accept.int) == 1 & is.null(antilog))
{
x.polygon=c(min(x),min(x),max(x),max(x))
y.polygon=c(diff.H0-accept.int/100*min(x),diff.H0+accept.int/100*min(x),diff.H0+accept.int/100*max(x),diff.H0-accept.int/100*max(x))
}
if(accept.int.perc == TRUE & length(accept.int) == 2 & is.null(antilog))
{
x.polygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
y.polygon=c(diff.H0-accept.int[1]/100*min(x),diff.H0+accept.int[1]/100*min(x),diff.H0+accept.int[1]/100*accept.int.threshold,diff.H0+accept.int[2]/100*accept.int.threshold,diff.H0+accept.int[2]/100*max(x),diff.H0-accept.int[2]/100*max(x),diff.H0-accept.int[2]/100*accept.int.threshold,diff.H0-accept.int[1]/100*accept.int.threshold)
}

	data=c(y,y.polygon)
	if("CI" %in% graph.int) data=c(data,Pred.CBLS[,CI1],Pred.CBLS[,CI2])
	if("CB" %in% graph.int) data=c(data,Pred.CBLS[,CB1],Pred.CBLS[,CB2])
	if("PI" %in% graph.int) data=c(data,Pred.CBLS[,PI1],Pred.CBLS[,PI2])
	if("GI" %in% graph.int) data=c(data,Pred.CBLS[,GI1],Pred.CBLS[,GI2])

	if(!("xlim" %in% names.arguments))
	{
	min_x=min(x)
	max_x=max(x)
	} else {
	min_x=arguments$xlim[1]
	max_x=arguments$xlim[2]
	}
	if(!("ylim" %in% names.arguments))
	{
	if(include.int == FALSE)
	{
	min_y=min(y)
	max_y=max(y)
	} else {
	min_y=min(data,na.rm=TRUE)
	max_y=max(data,na.rm=TRUE)
	}
	leg=(max_y-min_y)*0.05
	min_y=min_y-leg
	max_y=max_y+leg
	} else {
	min_y=arguments$ylim[1]
	max_y=arguments$ylim[2]
	}
	if(!("xlab" %in% names.arguments)) xname=x.title else xname=arguments$xlab
	if(!("ylab" %in% names.arguments)) yname=y.title else yname=arguments$ylab
	arguments.new=arguments
	arguments.new[c("xlim","ylim","xlab","ylab")]=NULL

	do.call(plot,args=c(list(x=1:10),list(y=1:10),list(type="n"),list(xlim=c(min_x,max_x)),list(ylim=c(min_y,max_y)),list(xlab=xname),list(ylab=yname),arguments.new))
	polygon(x.polygon,y.polygon,col="grey90",border=FALSE)
	
points(x,y,pch=19)
abline(a=diff.H0,b=0,lty=3)
lines(Pred.CBLS$X0,Pred.CBLS$Ypred,col=col.CBLS,lty=lty.CBLS)

cols=col.CBLS
ltys=lty.CBLS
texts="CBLS"

if("CI" %in% graph.int)
{
lines(Pred.CBLS$X0,Pred.CBLS[,CI1],col=col.CI,lty=lty.CI)
lines(Pred.CBLS$X0,Pred.CBLS[,CI2],col=col.CI,lty=lty.CI)
cols=c(cols,col.CI)
ltys=c(ltys,lty.CI)
texts=c(texts,paste0(conf.level,"% CI"))
}
if("CB" %in% graph.int)
{
lines(Pred.CBLS$X0,Pred.CBLS[,CB1],col=col.CB,lty=lty.CB)
lines(Pred.CBLS$X0,Pred.CBLS[,CB2],col=col.CB,lty=lty.CB)
cols=c(cols,col.CB)
ltys=c(ltys,lty.CB)
texts=c(texts,paste0(conf.level,"% CB"))
}
if("PI" %in% graph.int)
{
	if(sum(is.na(Pred.CBLS[,PI1])) == nrow(Pred.CBLS))
	{
	warning("No predictive intervals found to be displayed on the plot")
	} else {
	lines(Pred.CBLS$X0,Pred.CBLS[,PI1],col=col.PI,lty=lty.PI)
	lines(Pred.CBLS$X0,Pred.CBLS[,PI2],col=col.PI,lty=lty.PI)
	cols=c(cols,col.PI)
	ltys=c(ltys,lty.PI)
	texts=c(texts,paste0(pred.level,"% PI"))
	}
}
if("GI" %in% graph.int)
{
	if(sum(is.na(Pred.CBLS[,GI1])) == nrow(Pred.CBLS))
	{
	warning("No generalized predictive intervals found to be displayed on the plot")
	} else {
	lines(Pred.CBLS$X0,Pred.CBLS[,GI1],col=col.GI,lty=lty.GI)
	lines(Pred.CBLS$X0,Pred.CBLS[,GI2],col=col.GI,lty=lty.GI)
	cols=c(cols,col.GI)
	ltys=c(ltys,lty.GI)
	texts=c(texts,paste0(pred.level,"% GI"))
	}
}

if(sum(accept.int != 0) > 0) legend("topleft",c("Acceptance interval",leg0),bty="n",fill=c("grey90",0),lty=c(0,3),border=FALSE,horiz=TRUE) else legend("topleft",leg0,bty="n",lty=3,border=FALSE,horiz=TRUE)
legend("bottomright",lty=ltys,texts,col=cols,text.col=cols,bty="n",horiz=TRUE)
}#graph


} else {

res=results
if(!is.null(antilog)) res=antilog.pred(results=res,base=antilog)
x=res$data.MD$Means
y=res$data.MD$Differences
if(!is.null(accept.int.threshold)) {if(accept.int.threshold < min(x) | accept.int.threshold > max(x)) stop("accept.int.threshold should be between the minimum and maximum of (X+Y)/2 values, or use only one acceptance interval value.")}
Dmean=res$Intervals$'Mean differences'[1]
AI1=res$Intervals$Lower[1]
AI2=res$Intervals$Upper[1]
bTI1=res$Intervals$Lower[2]
bTI2=res$Intervals$Upper[2]
bgTI1=res$Intervals$Lower[3]
bgTI2=res$Intervals$Upper[3]

pred.level=res$Intervals$Meaning[1]
pred.level=substring(pred.level,1,regexpr("%",pred.level)[[1]]-1)
TI.conf.level=res$Intervals$Meaning[3]
TI.conf.level=substring(TI.conf.level,regexpr("h",TI.conf.level)[[1]]+2,regexpr("conf",TI.conf.level)[[1]]-3)

graph.MD=function()
{

if(!is.null(antilog)) diff.H0=1 else diff.H0=0

if(length(accept.int) == 1 & !is.null(antilog))
{
x.polygon=c(min(x),min(x),max(x),max(x))
y.polygon=c(diff.H0-accept.int/100,1/(diff.H0-accept.int/100),1/(diff.H0-accept.int/100),diff.H0-accept.int/100)
}
if(accept.int.perc == FALSE & length(accept.int) == 1 & is.null(antilog))
{
x.polygon=c(min(x),min(x),max(x),max(x))
y.polygon=c(diff.H0-accept.int,diff.H0+accept.int,diff.H0+accept.int,diff.H0-accept.int)
}
if(length(accept.int) == 2 & !is.null(antilog))
{
x.polygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
y.polygon=c(diff.H0-accept.int[1]/100,1/(diff.H0-accept.int[1]/100),1/(diff.H0-accept.int[1]/100),1/(diff.H0-accept.int[2]/100),1/(diff.H0-accept.int[2]/100),diff.H0-accept.int[2]/100,diff.H0-accept.int[2]/100,diff.H0-accept.int[1]/100)
}
if(accept.int.perc == FALSE & length(accept.int) == 2 & is.null(antilog))
{
x.polygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
y.polygon=c(diff.H0-accept.int[1],diff.H0+accept.int[1],diff.H0+accept.int[1],diff.H0+accept.int[2],diff.H0+accept.int[2],diff.H0-accept.int[2],diff.H0-accept.int[2],diff.H0-accept.int[1])
}
if(accept.int.perc == TRUE & length(accept.int) == 1 & is.null(antilog))
{
x.polygon=c(min(x),min(x),max(x),max(x))
y.polygon=c(diff.H0-accept.int/100*min(x),diff.H0+accept.int/100*min(x),diff.H0+accept.int/100*max(x),diff.H0-accept.int/100*max(x))
}
if(accept.int.perc == TRUE & length(accept.int) == 2 & is.null(antilog))
{
x.polygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
y.polygon=c(diff.H0-accept.int[1]/100*min(x),diff.H0+accept.int[1]/100*min(x),diff.H0+accept.int[1]/100*accept.int.threshold,diff.H0+accept.int[2]/100*accept.int.threshold,diff.H0+accept.int[2]/100*max(x),diff.H0-accept.int[2]/100*max(x),diff.H0-accept.int[2]/100*accept.int.threshold,diff.H0-accept.int[1]/100*accept.int.threshold)
}

	data=c(y,y.polygon)
	if("bTI" %in% graph.horiz.int) data=c(data,bTI1,bTI2)
	if("bgTI" %in% graph.horiz.int) data=c(data,bgTI1,bgTI2)

	if(!("xlim" %in% names.arguments))
	{
	min_x=min(x)
	max_x=max(x)
	} else {
	min_x=arguments$xlim[1]
	max_x=arguments$xlim[2]
	}
	if(!("ylim" %in% names.arguments))
	{
	if(include.int == FALSE)
	{
	min_y=min(y)
	max_y=max(y)
	} else {
	min_y=min(data)
	max_y=max(data)
	}
	leg=(max_y-min_y)*0.05
	min_y=min_y-leg
	max_y=max_y+leg
	} else {
	min_y=arguments$ylim[1]
	max_y=arguments$ylim[2]
	}
	if(!("xlab" %in% names.arguments)) xname=x.title else xname=arguments$xlab
	if(!("ylab" %in% names.arguments)) yname=y.title else yname=arguments$ylab
	arguments.new=arguments
	arguments.new[c("xlim","ylim","xlab","ylab")]=NULL

	do.call(plot,args=c(list(x=1:10),list(y=1:10),list(type="n"),list(xlim=c(min_x,max_x)),list(ylim=c(min_y,max_y)),list(xlab=xname),list(ylab=yname),arguments.new))
	polygon(x.polygon,y.polygon,col="grey90",border=FALSE)
	
points(x,y,pch=19)
abline(h=diff.H0,lty=3)
abline(h=Dmean)

if(sum(accept.int != 0) > 0) legend("topleft",c("Acceptance interval"),bty="n",fill="grey90",border=FALSE,horiz=TRUE)

if(is.null(antilog)) texte=c(expression(paste(bar(D))),leg0) else texte=c(expression(paste(bar(R))),leg0)
legend("topright",texte,bty="n",lty=c(1,3),border=FALSE,horiz=TRUE)

texte=NULL
ltys=NULL
cols=NULL

if("bTI" %in% graph.horiz.int)
{
abline(h=bTI1,col=col.bTI,lty=lty.bTI)
abline(h=bTI2,col=col.bTI,lty=lty.bTI)
texte=c(texte,as.expression(bquote(.(pred.level)*"%" ~ beta*TI)))
ltys=c(ltys,lty.bTI)
cols=c(cols,col.bTI)
}
if("bgTI" %in% graph.horiz.int)
{
abline(h=bgTI1,col=col.bgTI,lty=lty.bgTI)
abline(h=bgTI2,col=col.bgTI,lty=lty.bgTI)
texte=c(texte,as.expression(bquote(.(pred.level)*"%" ~ beta*gamma*TI ~ .(TI.conf.level)*"%")))
ltys=c(ltys,lty.bgTI)
cols=c(cols,col.bgTI)
}

if(!is.null(texte)) legend("bottomright",texte,lty=ltys,col=cols,text.col=cols,bty="n",horiz=TRUE)

}#graph

}#if else

graph.MD()

}
