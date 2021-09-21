#' @importFrom geiger treedata
fitDiscrete_return_lik <- function(phy, dat, k, model=c("ER","SYM","ARD","meristic"), ncores=1) {
    ## Function extracted from geiger. Returns only the likelihood for the model.
    td <- treedata(phy, dat)
    ## add check to make sure only unique data used
    if (nrow(td$data) != length(unique(rownames(td$data))))
        stop("Multiple records per tip label")

    phy=td$phy
    dat=td$data
    dd=dim(dat)
    trts=dd[2]
    ndat <- dat[,1]; # ah, gotta love R scoping...
    ## charStates <- sort(unique(ndat))
    charStates <- 1:k ## trying to allow for a non-strict model.
    constrain=model

    ## lik = geiger::mkn.lik(phy, ndat, constrain=constrain, transform="none")
    lik = mkn_lik_ext(phy = phy, dat = ndat, k = k, strict = FALSE, constrain=constrain, transform="none")
    attr(lik, "transform") = "none"
    argn=unlist(argn(lik))

    # translation of original character states
    attr(lik, "levels")<-charStates
    return( lik )
}
