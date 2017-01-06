#' @export
antilog.pred <-
function(results=NULL,base=10)
{
if(!(class(results) %in% c("BLS","CBLS","MD.horiz.lines","CIs.XY","CIs.MD"))) stop("results must be a BLS, CBLS, MD.horiz.lines, CIs.XY or CIs.MD class object (obtained after the BLS, CBLS, MD.horiz.lines function, FullCIs.XI or FullCIs.MD)")
if(base != 10 & base != "e") stop('base must be equal to "e" or 10')

if(base == 10) base = 10 else base = exp(1)

res=results

if(class(results) == "BLS")
{
res$Xij=base^(res$Xij)
res$Yik=base^(res$Yik)
res$Xi=base^(res$Xi)
res$Yi=base^(res$Yi)
res$Pred.BLS=as.data.frame(base^(res$Pred.BLS))
if(nrow(res$xpred.BLS) > 0) res$xpred.BLS=as.data.frame(base^(res$xpred.BLS))
}
if(class(results) == "CBLS")
{
res$Xij=base^(res$Xij)
res$Yik=base^(res$Yik)
res$Xi=base^(res$Xi)
res$Yi=base^(res$Yi)
res$Pred.CBLS=as.data.frame(base^(res$Pred.CBLS))
if(nrow(res$xpred.CBLS) > 0) res$xpred.CBLS=as.data.frame(base^(res$xpred.CBLS))
res$Mi=base^(res$Mi)
res$Di=base^(res$Di)
}
if(class(results) == "MD.horiz.lines")
{
res$data.MD=as.data.frame(base^(res$data.MD))
res$Intervals.horiz.lines[,1:3]=base^res$Intervals.horiz.lines[,1:3]
}
if(class(results) %in% c("CIs.XY","CIs.MD"))
{
res$Data.means=as.data.frame(base^res$Data.means)
res$Hyperbolic.intervals=base^res$Hyperbolic.intervals
}

return(res)
}
