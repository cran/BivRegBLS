#' @importFrom graphics abline legend lines par plot points polygon segments text title
#' @export
XY.plot <-
function(BLS.results=NULL,xname="X",yname="Y",antilog=NULL,accept.int=0,accept.int.perc=FALSE,accept.int.threshold=NULL,graph.int=c("CB","PI"),include.int=FALSE,col.BLS=1,col.CI=2,col.CB=3,col.PI=4,col.GI=5,lty.BLS=1,lty.CI=1,lty.CB=1,lty.PI=1,lty.GI=1,...)
{

if(!is.null(antilog)) {if(antilog != 10 & antilog != "e") stop('antilog must be equal to "e" or 10')}
if(class(BLS.results) != "BLS" & class(BLS.results) != "BLS.ht") stop("BLS.results must be a BLS or BLS.ht class object (obtained after the BLS or BLS.ht function)")
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

res=BLS.results
if(!is.null(antilog)) res=antilog.pred(results=res,base=antilog)
x=res$Xi
y=res$Yi
Pred.BLS=res$Pred.BLS
if(!is.null(accept.int.threshold)) {if(accept.int.threshold < min(x) | accept.int.threshold > max(x)) stop("accept.int.threshold should be between the minimum and maximum of X values, or use only one acceptance interval value.")}

coltitle=colnames(Pred.BLS)
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

arguments <- list(...)
names.arguments=names(arguments)

graph.XY=function()
{

	if (accept.int.perc == FALSE)
	{
		if(length(accept.int) == 1)
		{
		xpolygon=c(min(x),min(x),max(x),max(x))
		ypolygon=c(min(x)-accept.int,min(x)+accept.int,max(x)+accept.int,max(x)-accept.int)
		} else {
		xpolygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
		ypolygon=c(min(x)-accept.int[1],min(x)+accept.int[1],accept.int.threshold+accept.int[1],accept.int.threshold+accept.int[2],max(x)+accept.int[2],max(x)-accept.int[2],accept.int.threshold-accept.int[2],accept.int.threshold-accept.int[1])
		}
	} else {
		if(length(accept.int) == 1)
		{
		xpolygon=c(min(x),min(x),max(x),max(x))
		ypolygon=c(min(x)-accept.int/100*min(x),min(x)+accept.int/100*min(x),max(x)+accept.int/100*max(x),max(x)-accept.int/100*max(x))
		} else {
		xpolygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
		ypolygon=c(min(x)-accept.int[1]/100*min(x),min(x)+accept.int[1]/100*min(x),accept.int.threshold+accept.int[1]/100*accept.int.threshold,accept.int.threshold+accept.int[2]/100*accept.int.threshold,max(x)+accept.int[2]/100*max(x),max(x)-accept.int[2]/100*max(x),accept.int.threshold-accept.int[2]/100*accept.int.threshold,accept.int.threshold-accept.int[1]/100*accept.int.threshold)
		}
	}

	data=c(x,y,xpolygon,ypolygon)
	if("CI" %in% graph.int) data=c(data,Pred.BLS[,CI1],Pred.BLS[,CI2])
	if("CB" %in% graph.int) data=c(data,Pred.BLS[,CB1],Pred.BLS[,CB2])
	if("PI" %in% graph.int) data=c(data,Pred.BLS[,PI1],Pred.BLS[,PI2])
	if("GI" %in% graph.int) data=c(data,Pred.BLS[,GI1],Pred.BLS[,GI2])
	
	if(!("xlim" %in% names.arguments))
	{
	min_x=min(x,y)
	max_x=max(x,y)
	} else {
	min_x=arguments$xlim[1]
	max_x=arguments$xlim[2]
	}
	if(!("ylim" %in% names.arguments))
	{
	if(include.int == FALSE)
	{
	min_y=min_x
	max_y=max_x
	} else {
	min_y=min(data,na.rm=TRUE)
	max_y=max(data,na.rm=TRUE)
	}
	} else {
	min_y=arguments$ylim[1]
	max_y=arguments$ylim[2]
	}
	if(!("xlab" %in% names.arguments)) xname=xname else xname=arguments$xlab
	if(!("ylab" %in% names.arguments)) yname=yname else yname=arguments$ylab
	arguments.new=arguments
	arguments.new[c("xlim","ylim","xlab","ylab")]=NULL

	do.call(plot,args=c(list(x=1:10),list(y=1:10),list(type="n"),list(xlim=c(min_x,max_x)),list(ylim=c(min_y,max_y)),list(xlab=xname),list(ylab=yname),arguments.new))

	polygon(xpolygon,ypolygon,col="grey90",border=FALSE)

points(x,y,pch=19)
abline(a=0,b=1,lty=3)
lines(Pred.BLS$X0,Pred.BLS$Ypred,col=col.BLS,lty=lty.BLS)

cols=c(1,col.BLS)
ltys=c(3,lty.BLS)
texts=c("Y=X","BLS")

if("CI" %in% graph.int)
{
lines(Pred.BLS$X0,Pred.BLS[,CI1],col=col.CI,lty=lty.CI)
lines(Pred.BLS$X0,Pred.BLS[,CI2],col=col.CI,lty=lty.CI)
cols=c(cols,col.CI)
ltys=c(ltys,lty.CI)
texts=c(texts,paste0(conf.level,"% CI"))
}
if("CB" %in% graph.int)
{
lines(Pred.BLS$X0,Pred.BLS[,CB1],col=col.CB,lty=lty.CB)
lines(Pred.BLS$X0,Pred.BLS[,CB2],col=col.CB,lty=lty.CB)
cols=c(cols,col.CB)
ltys=c(ltys,lty.CB)
texts=c(texts,paste0(conf.level,"% CB"))
}
if("PI" %in% graph.int)
{
	if(sum(is.na(Pred.BLS[,PI1])) == nrow(Pred.BLS))
	{
	warning("No predictive intervals found to be displayed on the plot")
	} else {
	lines(Pred.BLS$X0,Pred.BLS[,PI1],col=col.PI,lty=lty.PI)
	lines(Pred.BLS$X0,Pred.BLS[,PI2],col=col.PI,lty=lty.PI)
	cols=c(cols,col.PI)
	ltys=c(ltys,lty.PI)
	texts=c(texts,paste0(pred.level,"% PI"))
	}
}
if("GI" %in% graph.int)
{
	if(sum(is.na(Pred.BLS[,GI1])) == nrow(Pred.BLS))
	{
	warning("No generalized predictive intervals found to be displayed on the plot")
	} else {
	lines(Pred.BLS$X0,Pred.BLS[,GI1],col=col.GI,lty=lty.GI)
	lines(Pred.BLS$X0,Pred.BLS[,GI2],col=col.GI,lty=lty.GI)
	cols=c(cols,col.GI)
	ltys=c(ltys,lty.GI)
	texts=c(texts,paste0(pred.level,"% GI"))
	}
}

if(sum(accept.int != 0) > 0) legend("topleft",c("Acceptance interval"),bty="n",fill="grey90",border=FALSE)
legend("bottomright",lty=ltys,texts,col=cols,text.col=cols,bty="n")
}

graph.XY()
}
