#' @importFrom grDevices dev.new
#' @export
GraphFullCIs.XY <-
function(FullCIs=NULL,BLS.estimate=NULL,lambda=NULL,xname="X",yname="Y",antilog=NULL,graph="joint.ellipse",accept.int=0,accept.int.perc=FALSE,accept.int.threshold=NULL,include.H0=TRUE,include.int=TRUE)
{

if(class(FullCIs) != "CIs.XY") stop("FullCIs must be an CIs.XY object obtained after the function FullCIs.XY")
if(!is.null(antilog)) {if(antilog != 10 & antilog != "e") stop('antilog must be equal to "e" or 10')}
if(!is.null(BLS.estimate) & class(BLS.estimate) != "BLS") stop("BLS.estimate must be a BLS class object (obtained after the BLS function)")
if(!is.null(lambda) & class(lambda) != "lambdas") stop("lambda must be a lambdas class object (obtained after the lambdas function)")
if(length(xname) > 1 | !is.character(xname)) stop("xname must be a character name")
if(length(yname) > 1 | !is.character(yname)) stop("yname must be a character name")
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
if(!is.null(accept.int.threshold)) if(length(accept.int.threshold) > 1 | !is.numeric(accept.int.threshold)) stop("accept.int.threshold must be one numeric value")
if(length(accept.int) == 2 & is.null(accept.int.threshold)) stop("If 2 values are used in the acceptance interval, then a threshold in X must be given by accept.int.threshold")

res=FullCIs
if(!is.null(antilog))
{
res=antilog.pred(results=res,base=antilog)
if(!is.null(BLS.estimate)) BLS.estimate=antilog.pred(results=BLS.estimate,base=antilog)
}

Ellipses.CB=res$Ellipses.CB
Slopes=res$Slopes
slopeNA=sum(is.na(Slopes[,3:4]))
if(slopeNA == 0) Slopesapprox=Slopes[,1:6] else Slopesapprox=Slopes[,c(1,2,5,6)]
Intercepts=res$Intercepts
ll=nrow(Slopes) # ll = number of lambda values
Hyperbolic.intervals=res$Hyperbolic.intervals
xx=Hyperbolic.intervals[,"X0",1]

lambdaXY1=lambda[[1]][2,3]
lambdaXY2=lambda[[1]][2,4]


graph_slopes=function()
{
	min_x=0.02
	max_x=50
	if(include.H0==TRUE) min_y=min(Slopesapprox[,-1],1,na.rm=TRUE) else min_y=min(Slopesapprox[,-1],na.rm=TRUE)
	if(include.H0==TRUE) max_y=max(Slopesapprox[,-1],1,na.rm=TRUE) else max_y=max(Slopesapprox[,-1],na.rm=TRUE)

	if(is.null(lambda)) kleg=0.02 else kleg=0.035
	leg=(max_y-min_y)*kleg
	min_y=min_y-leg
	max_y=max_y+leg

	plot(1:10,type="n",xlim=c(min_x,max_x),ylim=c(min_y,max_y),log="x",xlab=expression(lambda[XY]),ylab="Slope",cex.axis=1)
	if(!is.null(lambda))
	{
		polygon(c(lambdaXY1,lambdaXY1,lambdaXY2,lambdaXY2),c(min(Slopesapprox[,-1]),max(Slopesapprox[,-1]),max(Slopesapprox[,-1]),min(Slopesapprox[,-1])),col="grey90",border=FALSE)
		segments(lambda[[1]][2,1],min(Slopesapprox[,-1]),lambda[[1]][2,1],max(Slopesapprox[,-1]))
		text(lambda[[1]][2,1],max(Slopesapprox[,-1]),expression(paste(hat(lambda)[XY])),pos=3)
		legend("topright", expression(paste("CI ",lambda[XY])),bty='n',fill="grey90",border=FALSE)
	}
	
	if(slopeNA == 0)
	{
	lines(Slopes[c(-1,-ll),"Ratio variances"],Slopes[c(-1,-ll),3],lty=2)
	lines(Slopes[c(-1,-ll),"Ratio variances"],Slopes[c(-1,-ll),4],lty=2)
	}
	lines(Slopes[c(-1,-ll),"Ratio variances"],Slopes[c(-1,-ll),5])
	lines(Slopes[c(-1,-ll),"Ratio variances"],Slopes[c(-1,-ll),6])

	if(slopeNA == 0)
	{
	points(rep(min_x,3),Slopes[1,2:4],pch="_")
	segments(min_x,as.numeric(Slopes[1,3]),min_x,as.numeric(Slopes[1,4]))
	} else {
	points(rep(min_x,3),Slopes[1,c(2,5,6)],pch="_")
	segments(min_x,as.numeric(Slopes[1,5]),min_x,as.numeric(Slopes[1,6]))
	}

	points(rep(max_x,3),Slopes[ll,2:4],pch="_")
	segments(max_x,as.numeric(Slopes[ll,3]),max_x,as.numeric(Slopes[ll,4]))

	abline(h=1,lwd=1,lty=3)
	text(min_x,as.numeric(Slopes[1,2]),"OLSh",pos=4)
	text(max_x,Slopes[ll,"Estimate"],"OLSv",pos=2)
	
	coltitle=colnames(Slopes)
	conf.level=substring(coltitle[grep("CI",coltitle)[1]],regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-2,regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-1)
	
	if(slopeNA == 0)
	{
	legend("bottomleft", c(paste0("DR (exact ",conf.level,"% CI)"),paste0("BLS (approx ",conf.level,"% CI)")), lty = c(2,1),bty='n',horiz=TRUE)
	} else {
	legend("bottomleft", c(paste0("BLS (approx ",conf.level,"% CI)")), lty = c(1),bty='n',horiz=TRUE)
	}
}

graph_intercepts=function()
{
	min_x=0.02
	max_x=50
	if(include.H0==TRUE) min_y=min(Intercepts[,2:4],0) else min_y=min(Intercepts[,2:4])
	if(include.H0==TRUE) max_y=max(Intercepts[,2:4],0) else max_y=max(Intercepts[,2:4])

	if(is.null(lambda)) kleg=0.02 else kleg=0.035
	leg=(max_y-min_y)*kleg
	min_y=min_y-leg
	max_y=max_y+leg

	plot(1:10,type="n",xlim=c(min_x,max_x),ylim=c(min_y,max_y),log="x",xlab=expression(lambda[XY]),ylab="Intercept",cex.axis=1)
	if(!is.null(lambda))
	{
		polygon(c(lambdaXY1,lambdaXY1,lambdaXY2,lambdaXY2),c(min(Intercepts[,2:4]),max(Intercepts[,2:4]),max(Intercepts[,2:4]),min(Intercepts[,2:4])),col="grey90",border=FALSE)
		segments(lambda[[1]][2,1],min(Intercepts[,2:4]),lambda[[1]][2,1],max(Intercepts[,2:4]))
		text(lambda[[1]][2,1],max(Intercepts[,2:4]),expression(paste(hat(lambda)[XY])),pos=3)
		legend("topleft", expression(paste("CI ",lambda[XY])),bty='n',fill="grey90",border=FALSE)
	}

	lines(Intercepts[c(-1,-ll),"Ratio variances"],Intercepts[c(-1,-ll),3])
	lines(Intercepts[c(-1,-ll),"Ratio variances"],Intercepts[c(-1,-ll),4])

	points(rep(min_x,3),Intercepts[1,2:4],pch="_")
	segments(min_x,as.numeric(Intercepts[1,3]),min_x,as.numeric(Intercepts[1,4]))

	points(rep(max_x,3),Intercepts[ll,2:4],pch="_")
	segments(max_x,as.numeric(Intercepts[ll,3]),max_x,as.numeric(Intercepts[ll,4]))

	abline(h=0,lwd=1,lty=3)
	text(min_x,as.numeric(Intercepts[1,2]),"OLSh",pos=4)
	text(max_x,Intercepts[ll,"Estimate"],"OLSv",pos=2)

	coltitle=colnames(Intercepts)
	conf.level=substring(coltitle[grep("CI",coltitle)[1]],regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-2,regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-1)
	legend("bottomright", paste0("BLS (approx ",conf.level,"% CI)"), lty = 1,bty='n')
}

graph_ellipses.CB=function()
{
	x=res[[1]]$X
	if(sum(accept.int != 0) > 0 & is.null(antilog))
	{
	minx=min(x)
	maxx=max(x)
		if(accept.int.perc == FALSE)
		{
		slope_max=((maxx+accept.int[1])-(minx-accept.int[1]))/(maxx-minx)
		slope_min=((maxx-accept.int[1])-(minx+accept.int[1]))/(maxx-minx)
		intercept_max=minx+accept.int[1]-slope_min*minx
		intercept_min=minx-accept.int[1]-slope_max*minx
		xpolygon1=c(slope_min,1,slope_max,1)
		ypolygon1=c(intercept_max,accept.int[1],intercept_min,-accept.int[1])
			if(length(accept.int) > 1)
			{
			slope_max=((maxx+accept.int[2])-(minx-accept.int[2]))/(maxx-minx)
			slope_min=((maxx-accept.int[2])-(minx+accept.int[2]))/(maxx-minx)
			intercept_max=minx+accept.int[2]-slope_min*minx
			intercept_min=minx-accept.int[2]-slope_max*minx
			xpolygon2=c(slope_min,1,slope_max,1)
			ypolygon2=c(intercept_max,accept.int[2],intercept_min,-accept.int[2])
			}
		} else {
		slope_max=((maxx+accept.int[1]/100*maxx)-(minx-accept.int[1]/100*minx))/(maxx-minx)
		slope_min=((maxx-accept.int[1]/100*maxx)-(minx+accept.int[1]/100*minx))/(maxx-minx)
		intercept_max=minx+accept.int[1]/100*minx-slope_min*minx
		intercept_min=minx-accept.int[1]/100*minx-slope_max*minx
		xpolygon1=c(slope_min,(maxx+accept.int[1]/100*maxx)/(maxx),slope_max,(minx-accept.int[1]/100*minx)/minx)
		ypolygon1=c(intercept_max,0,intercept_min,0)
			if(length(accept.int) > 1)
			{
			slope_max=((maxx+accept.int[2]/100*maxx)-(minx-accept.int[2]/100*minx))/(maxx-minx)
			slope_min=((maxx-accept.int[2]/100*maxx)-(minx+accept.int[2]/100*minx))/(maxx-minx)
			intercept_max=minx+accept.int[2]/100*minx-slope_min*minx
			intercept_min=minx-accept.int[2]/100*minx-slope_max*minx
			xpolygon2=c(slope_min,(maxx+accept.int[2]/100*maxx)/(maxx),slope_max,(minx-accept.int[2]/100*minx)/minx)
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
		min_x=min(Ellipses.CB[,1,],1)
		max_x=max(Ellipses.CB[,1,],1)
		min_y=min(Ellipses.CB[,2,],0)
		max_y=max(Ellipses.CB[,2,],0)
		}
	} else {
	min_x=min(Ellipses.CB[,1,])
	max_x=max(Ellipses.CB[,1,])
	min_y=min(Ellipses.CB[,2,])
	max_y=max(Ellipses.CB[,2,])
	}

	leg=(max_y-min_y)*0.02
	min_y=min_y-leg
	max_y=max_y+leg

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
	if(!is.null(BLS.estimate)) lines(BLS.estimate$Ellipse.BLS,col=2,lwd=2)
	abline(v=1,lty=3)
	abline(h=0,lty=3)
	text(min(Ellipses.CB[,1,]),max(Ellipses.CB[,2,]),"OLSv",pos=3)
	text(max(Ellipses.CB[,1,]),min(Ellipses.CB[,2,]),"OLSh",pos=1)
}

graph_CIs=function()
{
	x=res$Data.means$X
	y=res$Data.means$Y
	data=c(x,y)
	if(!is.null(accept.int.threshold)) {if(accept.int.threshold < min(x) | accept.int.threshold > max(x)) stop("accept.int.threshold should be between the minimum and maximum of X values, or use only one acceptance interval value.")}

	coltitle=colnames(Hyperbolic.intervals)
	conf.level1=substring(coltitle[grep("CI",coltitle)[1]],1,regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-1)
	CI1=paste(conf.level1,"% CI Lower",sep="")
	CI2=paste(conf.level1,"% CI Upper",sep="")

	hyp_CI_min=apply(Hyperbolic.intervals[,CI1,],1,min)
	hyp_CI_max=apply(Hyperbolic.intervals[,CI2,],1,max)

	if(include.int==TRUE) data=c(data,hyp_CI_min,hyp_CI_max)

	min_x=min(data)
	max_x=max(data)

	plot(1:10,type="n",xlim=c(min_x,max_x),ylim=c(min_x,max_x),xlab=xname,ylab=yname)

	if (accept.int.perc == FALSE)
	{
		if(length(accept.int) == 1)
		{
		xpolygon=c(min(x),min(x),max(x),max(x))
		ypolygon=c(min(x)-accept.int,min(x)+accept.int,max(x)+accept.int,max(x)-accept.int)
		polygon(xpolygon,ypolygon,col="grey90",border=FALSE)
		} else {
		xpolygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
		ypolygon=c(min(x)-accept.int[1],min(x)+accept.int[1],accept.int.threshold+accept.int[1],accept.int.threshold+accept.int[2],max(x)+accept.int[2],max(x)-accept.int[2],accept.int.threshold-accept.int[2],accept.int.threshold-accept.int[1])
		polygon(xpolygon,ypolygon,col="grey90",border=FALSE)
		}
	} else {
		if(length(accept.int) == 1)
		{
		xpolygon=c(min(x),min(x),max(x),max(x))
		ypolygon=c(min(x)-accept.int/100*min(x),min(x)+accept.int/100*min(x),max(x)+accept.int/100*max(x),max(x)-accept.int/100*max(x))
		polygon(xpolygon,ypolygon,col="grey90",border=FALSE)
		} else {
		xpolygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
		ypolygon=c(min(x)-accept.int[1]/100*min(x),min(x)+accept.int[1]/100*min(x),accept.int.threshold+accept.int[1]/100*accept.int.threshold,accept.int.threshold+accept.int[2]/100*accept.int.threshold,max(x)+accept.int[2]/100*max(x),max(x)-accept.int[2]/100*max(x),accept.int.threshold-accept.int[2]/100*accept.int.threshold,accept.int.threshold-accept.int[1]/100*accept.int.threshold)
		polygon(xpolygon,ypolygon,col="grey90",border=FALSE)
		}
	}
	points(x,y,pch=19)

	abline(a=0,b=1,lty=3)
	lines(xx,hyp_CI_min)
	lines(xx,hyp_CI_max)
	if(!is.null(BLS.estimate))
	{
		# Find the confidence level, to display on the plot
		coltitle=colnames(BLS.estimate$Pred.BLS)
		conf.level2=substring(coltitle[grep("CI",coltitle)[1]],1,regexpr("%",coltitle[grep("CI",coltitle)[1]])[[1]]-1)
		CI1=paste(conf.level2,"% CI Lower",sep="")
		CI2=paste(conf.level2,"% CI Upper",sep="")

		#abline(a=as.numeric(BLS.estimate$Estimate.BLS["Intercept","Estimate"]),b=as.numeric(BLS.estimate$Estimate.BLS["Slope","Estimate"]),col=2)
		lines(BLS.estimate$Pred.BLS$X0,BLS.estimate$Pred.BLS$Ypred,col=2)
		lines(BLS.estimate$Pred.BLS$X0,BLS.estimate$Pred.BLS[,CI1],col=2,lwd=1)
		lines(BLS.estimate$Pred.BLS$X0,BLS.estimate$Pred.BLS[,CI2],col=2,lwd=1)
		legend("bottomright",lty=c(3,1,1),c("Y=X",paste0("Extreme ",conf.level1,"% CI"),paste0("BLS and ",conf.level2,"% CI")),col=c(1,1,2),bty="n")
	} else {
		legend("bottomright",lty=c(3,1),c("Y=X",paste0("Extreme ",conf.level1,"% CI")),bty="n")
	}
	if(sum(accept.int != 0) > 0) legend("topleft",c("Acceptance interval"),bty="n",fill="grey90",border=FALSE)
}


graph_CBs=function()
{
	x=res$Data.means$X
	y=res$Data.means$Y
	data=c(x,y)
	if(!is.null(accept.int.threshold)) {if(accept.int.threshold < min(x) | accept.int.threshold > max(x)) stop("accept.int.threshold should be between the minimum and maximum of X values, or use only one acceptance interval value.")}
	
	coltitle=colnames(Hyperbolic.intervals)
	conf.level1=substring(coltitle[grep("CB",coltitle)[1]],1,regexpr("%",coltitle[grep("CB",coltitle)[1]])[[1]]-1)
	CB1=paste(conf.level1,"% CB Lower",sep="")
	CB2=paste(conf.level1,"% CB Upper",sep="")

	hyp_CB_min=apply(Hyperbolic.intervals[,CB1,],1,min)
	hyp_CB_max=apply(Hyperbolic.intervals[,CB2,],1,max)

	if(include.int==TRUE) data=c(data,hyp_CB_min,hyp_CB_max)

	min_x=min(data)
	max_x=max(data)

	plot(1:10,type="n",xlim=c(min_x,max_x),ylim=c(min_x,max_x),xlab=xname,ylab=yname)

	if (accept.int.perc == FALSE)
	{
		if(length(accept.int) == 1)
		{
		xpolygon=c(min(x),min(x),max(x),max(x))
		ypolygon=c(min(x)-accept.int,min(x)+accept.int,max(x)+accept.int,max(x)-accept.int)
		polygon(xpolygon,ypolygon,col="grey90",border=FALSE)
		} else {
		xpolygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
		ypolygon=c(min(x)-accept.int[1],min(x)+accept.int[1],accept.int.threshold+accept.int[1],accept.int.threshold+accept.int[2],max(x)+accept.int[2],max(x)-accept.int[2],accept.int.threshold-accept.int[2],accept.int.threshold-accept.int[1])
		polygon(xpolygon,ypolygon,col="grey90",border=FALSE)
		}
	} else {
		if(length(accept.int) == 1)
		{
		xpolygon=c(min(x),min(x),max(x),max(x))
		ypolygon=c(min(x)-accept.int/100*min(x),min(x)+accept.int/100*min(x),max(x)+accept.int/100*max(x),max(x)-accept.int/100*max(x))
		polygon(xpolygon,ypolygon,col="grey90",border=FALSE)
		} else {
		xpolygon=c(min(x),min(x),accept.int.threshold,accept.int.threshold,max(x),max(x),accept.int.threshold,accept.int.threshold)
		ypolygon=c(min(x)-accept.int[1]/100*min(x),min(x)+accept.int[1]/100*min(x),accept.int.threshold+accept.int[1]/100*accept.int.threshold,accept.int.threshold+accept.int[2]/100*accept.int.threshold,max(x)+accept.int[2]/100*max(x),max(x)-accept.int[2]/100*max(x),accept.int.threshold-accept.int[2]/100*accept.int.threshold,accept.int.threshold-accept.int[1]/100*accept.int.threshold)
		polygon(xpolygon,ypolygon,col="grey90",border=FALSE)
		}
	}

	points(x,y,pch=19)

	abline(a=0,b=1,lty=3)
	lines(xx,hyp_CB_min)
	lines(xx,hyp_CB_max)
	if(!is.null(BLS.estimate))
	{
		coltitle=colnames(BLS.estimate$Pred.BLS)
		conf.level2=substring(coltitle[grep("CB",coltitle)[1]],1,regexpr("%",coltitle[grep("CB",coltitle)[1]])[[1]]-1)
		CB1=paste(conf.level2,"% CB Lower",sep="")
		CB2=paste(conf.level2,"% CB Upper",sep="")

		lines(BLS.estimate$Pred.BLS$X0,BLS.estimate$Pred.BLS$Ypred,col=2)
		lines(BLS.estimate$Pred.BLS$X0,BLS.estimate$Pred.BLS[,CB1],col=2,lwd=1)
		lines(BLS.estimate$Pred.BLS$X0,BLS.estimate$Pred.BLS[,CB2],col=2,lwd=1)
		legend("bottomright",lty=c(3,1,1),c("Y=X",paste0("Extreme ",conf.level1,"% CB"),paste0("BLS and ",conf.level2,"% CB")),col=c(1,1,2),bty="n")
	} else {
		legend("bottomright",lty=c(3,1),c("Y=X",paste0("Extreme ",conf.level1,"% CB")),bty="n")
	}
	if(sum(accept.int != 0) > 0) legend("topleft",c("Acceptance interval"),bty="n",fill="grey90",border=FALSE)
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
