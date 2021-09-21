#' @importFrom geiger argn
mkn_lik_ext <- function(phy, dat, k, strict = FALSE, constrain=c("ER","SYM","ARD","meristic"), transform=c("none", "EB", "lambda", "kappa", "delta", "white"), ...){
    ## Function extracted from geiger to allow for a non_strict model to run.
    phy = reorder(x = phy, order = "postorder")

    # control object for make.mkn()
	ct = list(method="exp")
	if(ct$method!="exp") stop(paste("method",sQuote(ct$method),"is not currently supported",sep=" "))

    # primary cache
	# k<-nlevels(as.factor(dat))
    if(is.character(dat)) dat=structure(as.factor(dat), names=names(dat))
    if(is.factor(dat)){
        levels=levels(dat)
        dat=structure(as.integer(dat), names=names(dat))
    } else {
        levels=sort(unique(dat))
    }
    if(k==2) if(all(constrain=="SYM")) constrain="ER"
    control <- geiger:::.check.control.mkn(ct, k)
    cache <- geiger:::.make.cache.mkn(phy, dat, k, strict=strict, control=ct, ultrametric=FALSE)
	cache$ordering=attributes(cache$info$phy)$order

    # tree transforms
	trns=match.arg(transform, c("none", "EB","lambda", "kappa", "delta", "white"))

	FUN=switch(trns,
    none=geiger:::.null.cache(cache),
    EB=geiger:::.eb.cache(cache),
    lambda=geiger:::.lambda.cache(cache),
    kappa=geiger:::.kappa.cache(cache),
    delta=geiger:::.delta.cache(cache),
    white=white.mkn(cache$states))

	if(trns=="white") return(FUN)

	ll.mkn=function(cache, control, ...) {
		k <- cache$info$k
		f.pars <- geiger:::.make.pars.mkn(k)
		f.pij <- geiger:::.make.pij.mkn(cache$info, control)
		idx.tip <- cache$idx.tip
		n.tip <- cache$n.tip
		n <- length(cache$len)
		map <- t(sapply(1:k, function(i) (1:k) + (i - 1) * k))
		idx.tip <- cbind(c(map[cache$states, ]), rep(seq_len(n.tip), k))
		children.C <- geiger:::.toC.int(t(cache$children))
		order.C <- geiger:::.toC.int(cache$order)

		.ll.mkn.exp=function(q, pars, intermediates=FALSE, preset = NULL) { # based on diversitree:::make.all.branches.mkn.exp
			if(is.null(argn(FUN))) new=FUN() else new=FUN(q)

			len.uniq <- sort(unique(new$len))
			len.idx <- match(new$len, len.uniq)

			if (!is.null(preset)) stop("Preset values not allowed")
			pij <- f.pij(len.uniq, pars)[, len.idx]
			lq <- numeric(n)
			branch.init <- branch.base <- matrix(NA, k, n)
			storage.mode(branch.init) <- "numeric"
			ans <- matrix(pij[idx.tip], n.tip, k)
			q <- rowSums(ans)
			branch.base[, seq_len(n.tip)] <- t.default(ans/q)
			lq[seq_len(n.tip)] <- log(q)
			ans <- .C("r_mkn_core", k = as.integer(k), n = length(order.C) -
            1L, order = order.C, children = children.C, pij = pij,
            init = branch.init, base = branch.base, lq = lq,
            NAOK = TRUE, PACKAGE="geiger")


			list(init = ans$init, base = ans$base, lq = ans$lq, vals = ans$init[, cache$root], pij = pij)
		}

        # build likelihood function
		attb=c(argn(FUN), cache$info$argn)
        rt=function(root="obs", root.p=NULL){
                    return(list(root=root, root.p=root.p))
        }
		if(is.null(argn(FUN))){ # NO TRANSFORM
			ll=function(pars, ...){
                rx=rt(...)
				qmat=f.pars(pars)
				ans=.ll.mkn.exp(q=NULL, pars=qmat, intermediates=FALSE)
				geiger:::.rootfunc.mkn(ans, qmat, root=rx$root, root.p=rx$root.p, intermediates=FALSE)
			}
		} else {
			ll=function(pars, ...){ # TREE TRANSFORM
                rx=rt(...)
				qmat=f.pars(pars[-1])
				ans=.ll.mkn.exp(q=pars[1], pars=qmat, intermediates=FALSE)
				geiger:::.rootfunc.mkn(ans, qmat, root=rx$root, root.p=rx$root.p, intermediates=FALSE)
			}

		}
		class(ll) <- c("mkn", "dtlik", "function")
		attr(ll,"argn") <- attb
        if(!is.null(levels)) attr(ll, "levels")=levels
		return(ll)
	}

	tmp=ll.mkn(cache, control)

    ## CONSTRAINTS
	if(!all(constrain=="ARD")){
		if(is.character(constrain)){
			cc=match.arg(constrain, c("ER","SYM","ARD","meristic"))
			tmp=geiger:::.constrain.k(tmp, model=cc, ...)
		} else {
			if(is.matrix(constrain)){
			    ## if(ncol(constrain)==max(dat)){ ## This check is preventing the constrain to work.
				## Modified the check to search for the number of states.
			    if(ncol(constrain)==k){
					tmp=geiger:::.constrain.m(tmp, m=constrain)
				}
			} else {
				stop("'constrain' must be supplied as a dummy matrix representing constraints on transition classes")
			}
		}
	}
	lik=function(pars, ...){
        if(missing(pars)) stop(paste("The following 'pars' are expected:\n\t", paste(argn(tmp), collapse="\n\t", sep=""), sep=""))

		pars=geiger:::.repars(pars, argn(tmp))
		tmp(pars, ...)
	}
	attributes(lik)<-attributes(tmp)
    attr(lik, "trns")<-!argn(lik)%in%argn(FUN)
	lik
}
