#' @importFrom phytools Map.Overlap
calcD <- function(t1, t2) {
  ## This computes the quantities to get the D metric.
  ## Depends only on the realization of the stochastic map. Do not seems to depend on the model of evolution.
	Do  <-  mapply(Map.Overlap,t1,t2,SIMPLIFY=FALSE)
	foo <- function(M){
		m1 <- rowSums(M)
		m2 <- colSums(M)
		as.matrix(m1)%*%t(m2)
	}
	De <- lapply(Do,foo)
	Dij <- mapply("-",Do,De,SIMPLIFY=FALSE)
	D <- sapply(Dij,function(x) sum(abs(x)))
	E_DX <- mean(D)
	E_Dij <- Reduce("+",Dij)/length(t1)
	list(E_DX=E_DX,E_Dij=E_Dij)
}
