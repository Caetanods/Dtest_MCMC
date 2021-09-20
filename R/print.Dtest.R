print.Dtest <-
function(x,...){
	if(hasArg(digits)) digits <- list(...)$digits
    else digits <- 4
    x <- lapply(x,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
	cat("Summary of results from D-test:\n")
	cat(paste(" E(D|X) = ",x$'E(D|X)',", P(D) = ",x$'P(D)',"\n",sep=""))
	cat("\n(Type ...$'E(Dij)' and ...$'P(Dij)' for\n pairwise E(D) and P-values.)\n")
}
