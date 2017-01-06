#' @export
df.WS <-
function(variances=NULL,k=rep(1,length(variances)),dfs=NULL)
{
if(is.null(variances) | is.null(dfs)) stop("No variances or no degrees of freedom are given")
if(length(variances) != length(dfs) | length(variances) != length(k) | length(dfs) != length(k)) stop("Lengths of 'variances', 'dfs' (degrees of freedom) or 'k' are different")
if(length(variances) == 1) stop("Only one value of variance is given, need at least two")
if(class(variances) != "numeric" | sum(variances <= 0) > 0) stop("Variances must be numeric and positive")
if(class(dfs) != "numeric" | sum(dfs <= 0) > 0) stop("Degrees of freedom (dfs) must be numeric and positive")
if(class(k) != "numeric") stop("k must be numeric")

Numerator=(sum(k*variances))^2
Denominator=sum((k*variances)^2/dfs)
df=Numerator/Denominator
return(df)
}
