#' @export
GraphFullCIs.MD <-
function(FullCIs=NULL,CBLS.estimate=NULL,lambda=NULL,xname="X",yname="Y",antilog=NULL,graph="joint.ellipse",accept.int=0,accept.int.perc=FALSE,accept.int.threshold=NULL,include.H0=TRUE,include.int=TRUE)
{

if(class(FullCIs) != "CIs.MD") stop("FullCIs must be an CIs.MD object obtained after the function FullCIs.MD")
if(!is.null(antilog)) warning("accept.int.perc is ignored if an antilog is applied")
if(!is.null(CBLS.estimate) & class(CBLS.estimate) != "CBLS") stop("CBLS.estimate must be a CBLS class object (obtained after the CBLS function)")
if(!is.null(lambda) & class(lambda) != "lambdas") stop("lambda must be a lambdas class object (obtained after the lambdas function)")
if(length(xname) > 1 | !is.character(xname)) stop("xname must be a character name")
if(length(yname) > 1 | !is.character(yname)) stop("yname must be a character name")
if(!is.null(antilog)) {if(antilog != 10 & antilog != "e") stop('antilog must be equal to "e" or 10')}
if(length(graph) > 1) stop("Only one type of graph is allowed in the argument 'graph'")
if(!(graph %in% c("slope","intercept","joint.ellipse","CI","CB","all"))) stop('graph must be equal at least to one of these values: "slope", "intercept", "joint.ellipse", "CI", "CB", "all"')

if(length(accept.int) > 2) stop("maximum two values are allowed in accept.int")
if(!is.numeric(accept.int)) stop("accept.int must be numeric")
if(length(accept.int.perc) > 1) stop("accept.int.perc must be TRUE or FALSE")
if(!(accept.int.perc %in% c(TRUE,FALSE))) stop("accept.int.perc must be TRUE or FALSE")
if(length(include.H0) > 1) stop("include.H0 must be TRUE or FALSE")
if(!(include.H0 %in% c(TRUE,FALSE))) stop("include.H0 must be TRUE or FALSE")
if(length(include.int) > 1) stop("include.int must be TRUE or FALSE")
if(!(include.int %in% c(TRUE,FALSE))) stop("include.int must be TRUE or FALSE")
if(!is.null(accept.int.threshold)) {if(length(accept.int.threshold) > 1 | !is.numeric(accept.int.threshold)) stop("accept.int.threshold must be one numeric value")}
if(length(accept.int) == 2 & is.null(accept.int.threshold)) stop("If 2 values are used in the acceptance interval, then a threshold in X must be given by accept.int.threshold")

res=FullCIs
if(!is.null(antilog))
{
res=antilog.pred(results=res,base=antilog)
if(!is.null(CBLS.estimate)) CBLS.estimate=antilog.pred(results=CBLS.estimate,base=antilog)
}

Ellipses.CB=res$Ellipses.CB
Slopes=res$Slopes
Intercepts=res$Intercepts
ll=nrow(Slopes) # ll = number of lambda values
Hyperbolic.intervals=res$Hyperbolic.intervals
xx=Hyperbolic.intervals[,"X0",1]

lambdaXY1=lambda[[1]][2,3]
lambdaXY2=lambda[[1]][2,4]
rho.hat=(lambda[[1]][2,1]-1)/(lambda[[1]][2,1]+1)
rho1=(lambdaXY1-1)/(lambdaXY1+1)
rho2=(lambdaXY2-1)/(lambdaXY2+1)

graph_slopes=function()
{
	min_x=-1
	max_x=1
	if(include.H0==TRUE) min_y=min(Slopes[,3:5],0) else min_y=min(Slopes[,3:5])
	if(include.H0==TRUE) max_y=max(Slopes[,3:5],0) else max_y=max(Slopes[,3:5])

	if(is.null(lambda)) kleg=0.02 else kleg=0.035
	leg=(max_y-min_y)*kleg
	min_y=min_y-leg
	max_y=max_y+leg

	plot(1:10,type="n",xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab=expression(rho[MD]),ylab="Slope",cex.axis=1)
	if(!is.null(lambda))
	{
		polygon(c(rho1,rho1,rho2,rho2),c(min(Slopes[,3:5]),max(Slopes[,3:5]),max(Slopes[,3:5]),min(Slopes[,3:5])),col="grey90",border=FALSE)
		segments(rho.hat,min(Slopes[,3:5]),rho.hat,max(Slopes[,3:5]))
		text(rho.hat,max(Slopes[,3:5]),expression(paste(hat(rho)[MD])),pos=3)
		legend("topright", expression(paste("CI ",rho[MD])),bty='n',fill="grey90",border=FALSE)
	}
	lines(Slopes[,"Correlation MD"],Slopes[,4],lwd=2)
	lines(Slopes[,"Correlation MD"],Slopes[,"Estimate"])
	lines(Slopes[,"Correlation MD"],Slopes[,5],lwd=2)
	abline(h=0,lwd=1,lty=3)

	coltitle=colnames(Slopes)
	conf.level=substring(coltitle[grep("CI",coltitle)[1]],regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-2,regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-1)
	legend("bottomleft", c("CBLS",paste0(conf.level,"% CI")), lty = c(1,1),lwd=c(1,2),bty='n',horiz=TRUE)
}

graph_intercepts=function()
{
	min_x=-1
	max_x=1
	if(include.H0==TRUE) min_y=min(Intercepts[,3:5],0) else min_y=min(Intercepts[,3:5])
	if(include.H0==TRUE) max_y=max(Intercepts[,3:5],0) else max_y=max(Intercepts[,3:5])

	if(is.null(lambda)) kleg=0.02 else kleg=0.035
	leg=(max_y-min_y)*kleg
	min_y=min_y-leg
	max_y=max_y+leg

	plot(1:10,type="n",xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab=expression(rho[MD]),ylab="Intercept",cex.axis=1)
	if(!is.null(lambda))
	{
		polygon(c(rho1,rho1,rho2,rho2),c(min(Intercepts[,3:5]),max(Intercepts[,3:5]),max(Intercepts[,3:5]),min(Intercepts[,3:5])),col="grey90",border=FALSE)
		segments(rho.hat,min(Intercepts[,3:5]),rho.hat,max(Intercepts[,3:5]))
		text(rho.hat,max(Intercepts[,3:5]),expression(paste(hat(rho)[MD])),pos=3)
		legend("topleft", expression(paste("CI ",rho[MD])),bty='n',fill="grey90",border=FALSE)
	}
	lines(Intercepts[,"Correlation MD"],Intercepts[,4],lwd=2)
	lines(Intercepts[,"Correlation MD"],Intercepts[,"Estimate"])
	lines(Intercepts[,"Correlation MD"],Intercepts[,5],lwd=2)
	abline(h=0,lwd=1,lty=3)

	coltitle=colnames(Intercepts)
	conf.level=substring(coltitle[grep("CI",coltitle)[1]],regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-2,regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-1)
	legend("bottomleft", c("CBLS",paste0(conf.level,"% CI")), lty = c(1,1),lwd=c(1,2),bty='n',horiz=TRUE)
}

graph_ellipses.CB=function()
{
x=res$Data.means$M
y=res$Data.means$D

	if(sum(accept.int != 0) > 0 & is.null(antilog))
	{
	minx=min(x)
	maxx=max(x)

		if(accept.int.perc == FALSE)
		{
		slope_max=((accept.int[1])-(-accept.int[1]))/(maxx-minx)
		slope_min=((-accept.int[1])-(accept.int[1]))/(maxx-minx)
		intercept_max=accept.int[1]-slope_min*minx
		intercept_min=-accept.int[1]-slope_max*minx
		xpolygon1=c(slope_min,0,slope_max,0)
		ypolygon1=c(intercept_max,accept.int[1],intercept_min,-accept.int[1])
			if(length(accept.int) > 1)
			{
			slope_max=((accept.int[2])-(-accept.int[2]))/(maxx-minx)
			slope_min=((-accept.int[2])-(accept.int[2]))/(maxx-minx)
			intercept_max=accept.int[2]-slope_min*minx
			intercept_min=-accept.int[2]-slope_max*minx
			xpolygon2=c(slope_min,0,slope_max,0)
			ypolygon2=c(intercept_max,accept.int[2],intercept_min,-accept.int[2])
			}
		} else {
		slope_max=((accept.int[1]/100*maxx)-(-accept.int[1]/100*minx))/(maxx-minx)
		slope_min=((-accept.int[1]/100*maxx)-(accept.int[1]/100*minx))/(maxx-minx)
		intercept_max=accept.int[1]/100*minx-slope_min*minx
		intercept_min=-accept.int[1]/100*minx-slope_max*minx
		xpolygon1=c(slope_min,(accept.int[1]/100*maxx)/(maxx),slope_max,(-accept.int[1]/100*minx)/minx)
		ypolygon1=c(intercept_max,0,intercept_min,0)
			if(length(accept.int) > 1)
			{
			slope_max=((accept.int[2]/100*maxx)-(-accept.int[2]/100*minx))/(maxx-minx)
			slope_min=((-accept.int[2]/100*maxx)-(accept.int[2]/100*minx))/(maxx-minx)
			intercept_max=accept.int[2]/100*minx-slope_min*minx
			intercept_min=-accept.int[2]/100*minx-slope_max*minx
			xpolygon2=c(slope_min,(accept.int[2]/100*maxx)/(maxx),slope_max,(-accept.int[2]/100*minx)/minx)
			ypolygon2=c(intercept_max,0,intercept_min,0)
			}
		}
	}

	if(include.H0==TRUE)
	{
		if(sum(accept.int != 0) > 0 & is.null(antilog))
		{
		if(length(accept.int) == 1) min_x=min(Ellipses.CB[,1,],xpolygon1) else min_x=min(Ellipses.CB[,1,],xpolygon1,xpolygon2)
		if(length(accept.int) == 1) max_x=max(Ellipses.CB[,1,],xpolygon1) else max_x=max(Ellipses.CB[,1,],xpolygon1,xpolygon2)
		if(length(accept.int) == 1) min_y=min(Ellipses.CB[,2,],ypolygon1) else min_y=min(Ellipses.CB[,2,],ypolygon1,ypolygon2)
		if(length(accept.int) == 1) max_y=max(Ellipses.CB[,2,],ypolygon1) else max_y=max(Ellipses.CB[,2,],ypolygon1,ypolygon2)
		} else {
		min_x=min(Ellipses.CB[,1,],0)
		max_x=max(Ellipses.CB[,1,],0)
		min_y=min(Ellipses.CB[,2,],0)
		max_y=max(Ellipses.CB[,2,],0)
		}
	} else {
	min_x=min(Ellipses.CB[,1,])
	max_x=max(Ellipses.CB[,1,])
	min_y=min(Ellipses.CB[,2,])
	max_y=max(Ellipses.CB[,2,])
	}

	leg=(max_y-min_y)*0.03
	min_y=min_y-leg
	max_y=max_y+leg
	leg=(max_x-min_x)*0.015
	min_x=min_x-leg
	max_x=max_x+leg

	plot(1:10,type="n",xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab='Slope',ylab='Intercept')
	if(sum(accept.int != 0) > 0 & is.null(antilog))
	{
		if(length(accept.int) == 1)
		{
		polygon(xpolygon1,ypolygon1,lty=3)
		legend("bottomleft", "Acceptance Interval",bty='n',lty=3,border=FALSE)
		} else {
		polygon(xpolygon1,ypolygon1,lty=3)
		polygon(xpolygon2,ypolygon2,lty=4)
		legend("bottomleft", c(paste("Acceptance Interval (",accept.int[1],")",sep=""),paste("Acceptance Interval (",accept.int[2],")",sep="")),bty='n',lty=c(3,4),border=FALSE)
		}	
	}
	for (j in 1:ll) lines(Ellipses.CB[,,j],col="grey")
	lines(Ellipses.CB[,,1],col=1,lwd=2)
	lines(Ellipses.CB[,,ll],col=1,lwd=2)
	if(!is.null(CBLS.estimate)) lines(CBLS.estimate$Ellipse.CBLS,col=2,lwd=2)
	abline(v=0,lty=3)
	abline(h=0,lty=3)
	text(min(Ellipses.CB[,1,]),max(Ellipses.CB[,2,]),expression(paste(rho[MD],"=1")),pos=3)
	text(max(Ellipses.CB[,1,]),min(Ellipses.CB[,2,]),expression(paste(rho[MD],"=-1")),pos=1)
	
}

if(is.null(antilog))
{
y.title=paste("Differences: D = ",yname," - ",xname,sep="")
x.title=paste("Averages: M = (",xname," + ",yname,")/2",sep="")
} else {
y.title=paste("Ratios: R = ",yname," / ",xname,sep="")
x.title=paste("Geom. Means: {",xname,", ",yname,"}",sep="")
}

if(!is.null(antilog)) diff.H0=1 else diff.H0=0

graph_CIs=function()
{
x=res$Data.means$M
y=res$Data.means$D
if(!is.null(accept.int.threshold)) {if(accept.int.threshold < min(x) | accept.int.threshold > max(x)) stop("accept.int.threshold should be between the minimum and maximum of (X+Y)/2 values, or use only one acceptance interval value.")}

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

	coltitle=colnames(Hyperbolic.intervals)
	conf.level1=substring(coltitle[grep("CI",coltitle)[1]],1,regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-1)
	CI1=paste(conf.level1,"% CI Lower",sep="")
	CI2=paste(conf.level1,"% CI Upper",sep="")

	hyp_CI_min=apply(Hyperbolic.intervals[,CI1,],1,min)
	hyp_CI_max=apply(Hyperbolic.intervals[,CI2,],1,max)

	data=c(y)
	if(include.H0 == TRUE) data=c(data,y.polygon)
	if(include.int == TRUE) data=c(data,hyp_CI_min,hyp_CI_max)

	min_x=min(x)
	max_x=max(x)
	min_y=min(data)
	max_y=max(data)

	leg=(max_y-min_y)*0.05
	min_y=min_y-leg
	max_y=max_y+leg

	plot(1:10,type="n",xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab=x.title,ylab=y.title)
	polygon(x.polygon,y.polygon,col="grey90",border=FALSE)
	points(x,y,pch=19)

	abline(a=diff.H0,b=0,lty=3)
	lines(xx,hyp_CI_min)
	lines(xx,hyp_CI_max)
	
	if(!is.null(CBLS.estimate))
	{
		coltitle=colnames(CBLS.estimate$Pred.CBLS)
		conf.level2=substring(coltitle[grep("CI",coltitle)[1]],1,regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-1)
		CI1=paste(conf.level2,"% CI Lower",sep="")
		CI2=paste(conf.level2,"% CI Upper",sep="")

		#abline(a=as.numeric(CBLS.estimate$Estimate.CBLS["Intercept","Estimate"]),b=as.numeric(CBLS.estimate$Estimate.CBLS["Slope","Estimate"]),col=2)
		lines(CBLS.estimate$Pred.CBLS$X0,CBLS.estimate$Pred.CBLS$Ypred,col=2)
		lines(CBLS.estimate$Pred.CBLS$X0,CBLS.estimate$Pred.CBLS[,CI1],col=2,lwd=1)
		lines(CBLS.estimate$Pred.CBLS$X0,CBLS.estimate$Pred.CBLS[,CI2],col=2,lwd=1)
		legend("bottomright",lty=c(1,1),c(paste0("Extreme ",conf.level1,"% CI"),paste0("CBLS and ",conf.level2,"% CI")),col=c(1,2),bty="n",horiz=TRUE)
	} else {
		legend("bottomright",lty=1,c(paste0("Extreme ",conf.level1,"% CI")),col=1,bty="n",horiz=TRUE)
	}
	if(is.null(antilog)) leg0="D=0" else leg0="R=1"
	if(sum(accept.int != 0) > 0) legend("topleft",c("Acceptance interval",leg0),bty="n",fill=c("grey90",0),lty=c(0,3),border=FALSE,horiz=TRUE) else legend("topleft",leg0,bty="n",lty=3,border=FALSE,horiz=TRUE)
}


graph_CBs=function()
{

x=res$Data.means$M
y=res$Data.means$D
if(!is.null(accept.int.threshold)) {if(accept.int.threshold < min(x) | accept.int.threshold > max(x)) stop("accept.int.threshold should be between the minimum and maximum of (X+Y)/2 values, or use only one acceptance interval value.")}

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

	coltitle=colnames(Hyperbolic.intervals)
	conf.level1=substring(coltitle[grep("CB",coltitle)[1]],1,regexpr("%",coltitle[grep("CB",coltitle)[1]])[[1]]-1)
	CB1=paste(conf.level1,"% CB Lower",sep="")
	CB2=paste(conf.level1,"% CB Upper",sep="")

	hyp_CB_min=apply(Hyperbolic.intervals[,CB1,],1,min)
	hyp_CB_max=apply(Hyperbolic.intervals[,CB2,],1,max)

	data=c(y)
	if(include.H0 == TRUE) data=c(data,y.polygon)
	if(include.int == TRUE) data=c(data,hyp_CB_min,hyp_CB_max)

	min_x=min(x)
	max_x=max(x)
	min_y=min(data)
	max_y=max(data)

	leg=(max_y-min_y)*0.05
	min_y=min_y-leg
	max_y=max_y+leg

	plot(1:10,type="n",xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab=x.title,ylab=y.title)
	polygon(x.polygon,y.polygon,col="grey90",border=FALSE)
	points(x,y,pch=19)

	abline(a=diff.H0,b=0,lty=3)
	lines(xx,hyp_CB_min)
	lines(xx,hyp_CB_max)

	if(!is.null(CBLS.estimate))
	{
		coltitle=colnames(CBLS.estimate$Pred.CBLS)
		conf.level2=substring(coltitle[grep("CB",coltitle)[1]],1,regexpr("%",coltitle[grep("CB",coltitle)[1]])[[1]]-1)
		CB1=paste(conf.level2,"% CB Lower",sep="")
		CB2=paste(conf.level2,"% CB Upper",sep="")

		#abline(a=as.numeric(CBLS.estimate$Estimate.CBLS["Intercept","Estimate"]),b=as.numeric(CBLS.estimate$Estimate.CBLS["Slope","Estimate"]),col=2)
		lines(CBLS.estimate$Pred.CBLS$X0,CBLS.estimate$Pred.CBLS$Ypred,col=2)
		lines(CBLS.estimate$Pred.CBLS$X0,CBLS.estimate$Pred.CBLS[,CB1],col=2,lwd=1)
		lines(CBLS.estimate$Pred.CBLS$X0,CBLS.estimate$Pred.CBLS[,CB2],col=2,lwd=1)
		legend("bottomright",lty=c(1,1),c(paste0("Extreme ",conf.level1,"% CB"),paste0("CBLS and ",conf.level2,"% CB")),col=c(1,2),bty="n",horiz=TRUE)
	} else {
		legend("bottomright",lty=1,c(paste0("Extreme ",conf.level1,"% CB")),col=1,bty="n",horiz=TRUE)
	}
	if(is.null(antilog)) leg0="D=0" else leg0="R=1"
	if(sum(accept.int != 0) > 0) legend("topleft",c("Acceptance interval",leg0),bty="n",fill=c("grey90",0),lty=c(0,3),border=FALSE,horiz=TRUE) else legend("topleft",leg0,bty="n",lty=3,border=FALSE,horiz=TRUE)
}

if("all" %in% graph)
{
dev.new();graph_slopes()
dev.new();graph_intercepts()
dev.new();graph_ellipses.CB()
dev.new();graph_CIs()
dev.new();graph_CBs()
} else {
if("slope" %in% graph) graph_slopes()
if("intercept" %in% graph) graph_intercepts()
if("joint.ellipse" %in% graph) graph_ellipses.CB()
if("CI" %in% graph) graph_CIs()
if("CB" %in% graph) graph_CBs()
}

}
