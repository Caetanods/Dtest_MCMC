#' Function to perform the MCMC posterior predictive simulations
#'
#' @param t1 a list of stochastic maps.
#' @param t2 a list of stochastic maps.
#' @param tag a name (pattern) used to name the files produced by the simulation.
#' @param dir_files the directory where the files will be save. You WILL need to know where this location is!
#' @param ngen number of generations for each MCMC chain.
#' @param nsim number of posterior predictive simulation replicates.
#' @param from_to index of the posterior predictive replicates to run. Can be used to structure an analysis over multiple processors and/or computers.
#' @param read_only if the function should only read the existing files instead of producing the simulations.
#' @param model_mat1 the model for the trait 1. Same format as in the "model" parameter for the function "ace" in the package "ace".
#' @param model_mat2 the model for the trait 2. Same format as in the "model" parameter for the function "ace" in the package "ace".
#' @param root_p1 the vector of root probabilities for trait 1.
#' @param root_p2 the vector of root probabilities for trait 2.
#' @param prior_exp the parameter for the exponential prior for the rates of transition of the Markov model fitted to the data.
#' @param scale_mcmc the scaling factor for the Metropolis step of the MCMC. See function "metrop" in package "mcmc".
#'
#' @return write a large number of ".rds" files to the directory. Please make sure that you have read and write permissions to the "dir_files" location! Also make sure to record the location of the "dir_files" directory.
#' @export
Dtest_MCMC_informed_root <- function(t1, t2, tag = "test", dir_files = "./partials", ngen = 5000, nsim = 1000, from_to = NULL, read_only = FALSE, model_mat1, model_mat2, root_p1, root_p2, prior_exp = 0.5, scale_mcmc = 1){

    ## This is similar to the previous function. But in this case the root state is informed by a root probability vector.
    ## If read_only = TRUE no new simulation will be done and the test will use only the replicates already available in the 'dir_files'.
    cat("This test uses posterior predictive simulations. Please check the documentation.\n")
	cat("Running. (This WILL take some time.)\n")

	## Creates a directory to store the partial results.
	dir.create( dir_files, showWarnings = FALSE )

	## Make a list of the files already existent in the directory.
	file_list <- list.files(path = dir_files, pattern = "*.rds", full.names = FALSE)

	## If using 'read_only' option, then we need the number of observed replicates.
	## The expected number is 'length(t1)', but it could be lower.
	## This need to be computed before the loop starts.
	## If read_only, then check if the files exist, if not then skip.
	nrep <- length(t1) ## The expected number of replicates.
	if( read_only ){
	    check_simmap_t1_files <- length( grep(x = file_list, pattern = paste0("simmap_t1_rep", tag)) )
	    check_mcmc_t1_files <- length( grep(x = file_list, pattern = paste0("mcmc_t1_rep", tag)) )
	    check_simmap_t2_files <- length( grep(x = file_list, pattern = paste0("simmap_t2_rep", tag)) )
	    check_mcmc_t2_files <- length( grep(x = file_list, pattern = paste0("mcmc_t2_rep", tag)) )
	    obs_nrep <- min(check_simmap_t1_files, check_mcmc_t1_files
	                    , check_simmap_t2_files, check_mcmc_t2_files)
	} else{
	    obs_nrep <- nrep
	}

	## Need to append the class "multiPhylo" to the list object.
	if( !inherits(t1,"multiPhylo") ){
	  class(t1) <- c(class(t1), "multiPhylo" )
	}
	if( !inherits(t2,"multiPhylo") ){
	  class(t2) <- c(class(t2), "multiPhylo" )
	}

	## Get the number of states on each of the vectors.
	levs1 <- sort(unique(as.vector(phytools::getStates(t1,"tips"))))
	k1 <- length(levs1)
	levs2 <- sort(unique(as.vector(phytools::getStates(t2,"tips"))))
	k2 <- length(levs2)

	## Get the number of parameters for each of the models:
	npar_t1 <- sum(unique(c(model_mat1)) > 0, na.rm = TRUE)
	npar_t2 <- sum(unique(c(model_mat2)) > 0, na.rm = TRUE)

	## Calculates the D metric given the observed stochastic maps.
	obj <- calcD(t1,t2) ## Note that we have a distribution of stochastic maps here.
	E_DX <- obj$E_DX
	E_Dij <- obj$E_Dij

	## posterior prediction (Here is where we do the simulations)
	PD <- 0
	## Matrix to store the results.
	Pdij <- matrix(0,k1,k2,dimnames=list(levs1,levs2))
	cat("|")

	## Enable the function to only do the replicates that we want to.
	## Because the function will save all the partials, this is an efficient way to complete the replicates.
	if( is.null(from_to) ){
	    for_loop <- 1:nrep
	} else{
	    if( !is.numeric(from_to) ) stop("from_to needs to be a numeric vector.")
	    if( max( from_to ) > nrep ) stop( "from_to is larger than the number of stochastic maps available.")
	    for_loop <- from_to
	}

	for(i in for_loop){ ## For each of the posterior predictive reps we need to do a full MCMC.

	    ## If read_only, then check if the files exist, if not then skip.
	    if( read_only ){
	        check_files <- c(  paste0("mcmc_t1_rep"  , tag, i, ".rds") %in% file_list
	                         , paste0("simmap_t1_rep", tag, i, ".rds") %in% file_list
	                         , paste0("mcmc_t2_rep"  , tag, i, ".rds") %in% file_list
	                         , paste0("simmap_t2_rep", tag, i, ".rds") %in% file_list )
	        if( !all(check_files) ) next
	    }

	    ## Need to fix the rownames of the Q matrix:
	    Q_t1 <- t1[[i]]$Q
	    k_t1 <- ncol( Q_t1 )
	    rownames( Q_t1 ) <- colnames( Q_t1 )
	    ## Simulation under the model to make data for a posterior predictive. (For t1 and t2)
		## x <- phytools::sim.Mk(t1[[i]], Q_t1, anc = sample(x = colnames(Q_t1), size = 1, prob = root_p1))
		## Here using rejection sampling to make sure we have the correct number of states.
		x <- reject_sample_MK(tree = t1[[i]], Q = Q_t1, root_pi = root_p1)
		x <- setNames(as.character(x), names(x))
		Q_t2 <- t2[[i]]$Q
		k_t2 <- ncol( Q_t2 )
		rownames( Q_t2 ) <- colnames( Q_t2 )
		## y <- phytools::sim.Mk(t2[[i]], Q_t2, anc = sample(x = colnames(Q_t2), size = 1, prob = root_p2))
		y <- reject_sample_MK(tree = t2[[i]], Q = Q_t2, root_pi = root_p2)
		y <- setNames(as.character(y), names(y))

		## Here is the bottleneck in terms of the time to run the analysis.
		## I am going to use "diversitree" to perform the MCMC and "ratematrix" to make the Stochastic maps.

		########################################
		## Block for the first tree.
		########################################
		x_num <- sapply(x, function(xxx) which(levs1 == xxx)) ## Need to transform the vector into numbers.
		## Construct the likelihood function.
		mkn_t1 <- fitDiscrete_return_lik(phy = t1[[i]], dat = x_num, k = k_t1, model = model_mat1, ncores = 1)
		post_fn_t1 <- function(x){
		    x <- exp(x) ## Search in log space.
		    if( any(x < 0.0) ) return( -Inf )
		    log_prior <- sum( dexp(x = x, rate = prior_exp, log = TRUE) )
		    log_lik <- mkn_t1(pars = x, root = "given", root.p = root_p1)
		    if( !is.finite(log_lik + log_prior) ) return( -Inf )
		    return( log_lik + log_prior )
		}

		## Before running the MCMC, check if the MCMC file is available.
		if( paste0("mcmc_t1_rep", tag, i, ".rds") %in% file_list ){
			## Read the file.
			mcmc_t1 <- readRDS( paste0(dir_files, "/mcmc_t1_rep", tag, i, ".rds") )
		} else{
			## Run the MCMC and save the rds file.
			## Get the MLE as starting state.
		    t1_tmp <- ratematrix::mergeSimmap(phy = t1[[i]], drop.regimes = TRUE)
		    invisible(capture.output( start_mle_guess_corHMM <- corHMM(phy = t1_tmp, data = data.frame(spp = names(x_num), trait = as.numeric(x_num)), rate.cat = 1, rate.mat = model_mat1, root.p = root_p1) ))
		    start_mle_guess <- sapply(1:max(model_mat1, na.rm = TRUE), function(x) start_mle_guess_corHMM$solution[which(model_mat1 == x)][1])
		    if( any(start_mle_guess < 1e-6) ){
		        start_mle_guess[start_mle_guess < 1e-6] <- 1e-6
		    }
		    ## The code below works, but the guess is not good. corHMM does a much better job.
		#     start_mle_guess <- optim(par = runif(n = npar_t1, min = min.start, max = max.start)
		# 	                         , method = "L-BFGS-B"
		# 	                         , lower = rep(x = unif_prior1[1], times = npar_t1)
		# 	                         , upper = rep(x = unif_prior1[2], times = npar_t1)
		# 	                         , fn = mkn_t1, root = "given", root.p = root_p1
		# 	                         , control = list(fnscale = -1))
			## Run the MCMC.
		    mcmc_t1 <- mcmc::metrop(obj = post_fn_t1, initial = log(start_mle_guess), nbatch = ngen
		                            , blen = 1, nspac = 1, scale = scale_mcmc)
		    mcmc_t1$batch <- exp(mcmc_t1$batch) ## Need to untransform the parameter space.
			## Save the partial result.
			## Need to remember that this posterior is in a distinct format.
			saveRDS(mcmc_t1, file = paste0(dir_files, "/mcmc_t1_rep", tag, i, ".rds") )
		}

		## Sample 1000 Q matrices from the MCMC in order to make the stochastic maps.
		ss <- sample(x = seq(from = round(ngen/2), to = ngen, by = 1), size = nsim, replace = FALSE)
		dt_t1 <- mcmc_t1$batch[ss,]
		## Transform the sample in a list of Q matrices.
		Q_mcmc_t1 <- list()
		for( w in 1:nsim ){
		    ## Need to reconstruct using the model matrix.
		    mm <- model_mat1
		    for( j in 1:npar_t1 ) mm[model_mat1 == j] <- dt_t1[w,j]
		    diag( mm ) <- 0.0 ## Make sure the diag is 0!
		    diag( mm ) <- rowSums( mm ) * -1
		    rownames(mm) <- colnames(mm) <- 1:k_t1
		    Q_mcmc_t1[[w]] <- mm
		}

		## Before doing the stochastic maps, check if these were already done.
		if( paste0("simmap_t1_rep", tag, i, ".rds") %in% file_list ){
			## Read in the file and skip the stochastic maps.
			T1 <- readRDS( paste0(dir_files, "/simmap_t1_rep", tag, i, ".rds") )
		} else{
			## Make the stochastic maps and save the rds file.
			T1 <- lapply(Q_mcmc_t1, function(xxx) ratematrix::fastSimmap(tree = t1[[i]], x = x
							, pi = root_p1, Q = xxx, silence = TRUE, max_nshifts = 200))
			## Check for problems in the stochastic map simulation:
			t1_check_class <- sapply(T1, function(x) "phylo" %in% class(x))
			if( any(!t1_check_class) ){
			    ## Some have problems. Need to solve before moving on.
			    problem_id <- which(!t1_check_class)
			    for( solve_id in problem_id ){
			        init_nshifts <- 2000
			        n_try_t1_simmap <- 0
			        while( !"phylo" %in% class(T1[[solve_id]]) ){
			            T1[[solve_id]] <- ratematrix::fastSimmap(tree = t1[[i]], x = y
			                                                     , pi = root_p1
			                                                     , Q = Q_mcmc_t1[[solve_id]]
			                                                     , silence = TRUE
			                                                     , max_nshifts = init_nshifts)
			            init_nshifts <- init_nshifts * 10
			            n_try_t1_simmap <- n_try_t1_simmap + 1
			            if( n_try_t1_simmap > 5 ) break
			        }
			        if( !"phylo" %in% class(T1[[solve_id]]) ){
			            ## Still not good, give up and re-sample:
			            ss <- sample(x = seq(from = round(ngen/2), to = ngen, by = 1), size = 1)
			            Q_mcmc_t1_tmp <- model_mat1
			            for( j in 1:npar_t1 ) Q_mcmc_t1_tmp[model_mat1 == j] <- mcmc_t1$batch[ss,j]
			            diag( Q_mcmc_t1_tmp ) <- 0.0 ## Make sure the diag is 0!
			            diag( Q_mcmc_t1_tmp ) <- rowSums( Q_mcmc_t1_tmp ) * -1
			            rownames(Q_mcmc_t1_tmp) <- colnames(Q_mcmc_t1_tmp) <- 1:k_t1

			            ## Another round of sampling!
			            init_nshifts <- 2000
			            n_try_t1_simmap <- 0
			            while( !"phylo" %in% class(T1[[solve_id]]) ){
			                T1[[solve_id]] <- ratematrix::fastSimmap(tree = t1[[i]], x = y
			                                                         , pi = root_p1
			                                                         , Q = Q_mcmc_t1_tmp
			                                                         , silence = TRUE
			                                                         , max_nshifts = init_nshifts)
			                init_nshifts <- init_nshifts * 10
			                n_try_t1_simmap <- n_try_t1_simmap + 1
			                if( n_try_t1_simmap > 5 ) break ## Give up!
			            }
			        }
			    }
			}
			## Save the T1 samples.
			saveRDS(T1, file = paste0(dir_files, "/simmap_t1_rep", tag, i, ".rds") )
		}

		########################################
		## Block for the second tree.
		########################################
		y_num <- sapply(y, function(yyy) which(levs2 == yyy)) ## Need to transform the vector into numbers.
		mkn_t2 <- fitDiscrete_return_lik(phy = t2[[i]], dat = y_num, k = k_t2, model = model_mat2, ncores = 1)
		post_fn_t2 <- function(x){
		    x <- exp(x)
		    if( any(x < 0.0) ) return( -Inf )
		    log_prior <- sum( dexp(x = x, rate = prior_exp, log = TRUE) )
		    log_lik <- mkn_t2(pars = x, root = "given", root.p = root_p2)
		    if( !is.finite(log_lik + log_prior) ) return( -Inf )
		    return( log_lik + log_prior )
		}

		## Before running the MCMC, check if the MCMC file is available.
		if( paste0("mcmc_t2_rep", tag, i, ".rds") %in% file_list ){
			## Read the file.
			mcmc_t2 <- readRDS( paste0(dir_files, "/mcmc_t2_rep", tag, i, ".rds") )
		} else{
		    ## Run the MCMC and save the rds file.
		    ## Get the MLE as starting state.
		    t2_tmp <- ratematrix::mergeSimmap(phy = t2[[i]], drop.regimes = TRUE)
		    invisible(capture.output( start_mle_guess_corHMM <- corHMM(phy = t2_tmp, data = data.frame(spp = names(y_num), trait = as.numeric(y_num)), rate.cat = 1, rate.mat = model_mat2, root.p = root_p2) ))
		    start_mle_guess <- sapply(1:max(model_mat2, na.rm = TRUE), function(x) start_mle_guess_corHMM$solution[which(model_mat2 == x)][1])
		    if( any(start_mle_guess < 1e-6) ){
		        start_mle_guess[start_mle_guess < 1e-6] <- 1e-6
		    }
		    # start_mle_guess <- optim(par = runif(n = npar_t2, min = min.start, max = max.start)
		    #                          , method = "L-BFGS-B"
		    #                          , lower = rep(x = unif_prior2[1], times = npar_t2)
		    #                          , upper = rep(x = unif_prior2[2], times = npar_t2)
		    #                          , fn = mkn_t2, root = "given", root.p = root_p2
		    #                          , control = list(fnscale = -1))
		    ## Run the MCMC.
		    mcmc_t2 <- mcmc::metrop(obj = post_fn_t2, initial = log(start_mle_guess), nbatch = ngen
		                            , blen = 1, nspac = 1, scale = scale_mcmc)
		    mcmc_t2$batch <- exp(mcmc_t2$batch) ## Need to unstrasform the MCMC samples
			## Save the partial result.
			saveRDS(mcmc_t2, file = paste0(dir_files, "/mcmc_t2_rep", tag, i, ".rds") )
		}

		ss <- sample(x = seq(from = round(ngen/2), to = ngen, by = 1), size = nsim, replace = FALSE)
		dt_t2 <- mcmc_t2$batch[ss,]
		## Transform the sample in a list of Q matrices.
		Q_mcmc_t2 <- list()
		for( w in 1:nsim ){
		    ## Need to reconstruct using the model matrix.
		    mm <- model_mat2
		    for( j in 1:npar_t2 ) mm[model_mat2 == j] <- dt_t2[w,j]
		    diag( mm ) <- 0.0 ## Make sure the diag is 0!
		    diag( mm ) <- rowSums( mm ) * -1
		    rownames(mm) <- colnames(mm) <- 1:k_t2
		    Q_mcmc_t2[[w]] <- mm
		}

		## Before doing the stochastic maps, check if these were already done.
		if( paste0("simmap_t2_rep", tag, i, ".rds") %in% file_list ){
			## Read in the file and skip the stochastic maps.
			T2 <- readRDS( paste0(dir_files, "/simmap_t2_rep", tag, i, ".rds") )
		} else{
		    ## Make the stochastic maps and save the rds file.
		    T2 <- lapply(Q_mcmc_t2, function(xxx) ratematrix::fastSimmap(tree = t2[[i]], x = x
		                                                                 , pi = root_p2, Q = xxx
		                                                                 , silence = TRUE
		                                                                 , max_nshifts = 200))
		    ## Check for problems in the stochastic map simulation:
		    t2_check_class <- sapply(T2, function(x) "phylo" %in% class(x))
		    if( any(!t2_check_class) ){
		        ## Some have problems. Need to solve before moving on.
		        problem_id <- which(!t2_check_class)
		        for( solve_id in problem_id ){
		            init_nshifts <- 2000
		            n_try_t2_simmap <- 0
		            while( !"phylo" %in% class(T2[[solve_id]]) ){
		                T2[[solve_id]] <- ratematrix::fastSimmap(tree = t2[[i]], x = y
		                                                         , pi = root_p2
		                                                         , Q = Q_mcmc_t2[[solve_id]]
		                                                         , silence = TRUE
		                                                         , max_nshifts = init_nshifts)
		                init_nshifts <- init_nshifts * 10
		                n_try_t2_simmap <- n_try_t2_simmap + 1
		                if( n_try_t2_simmap > 5 ) break
		            }
		            if( !"phylo" %in% class(T2[[solve_id]]) ){
		                ## Still not good, give up and re-sample:
		                ss <- sample(x = seq(from = round(ngen/2), to = ngen, by = 1), size = 1)
		                Q_mcmc_t2_tmp <- model_mat2
		                for( j in 1:npar_t2 ) Q_mcmc_t2_tmp[model_mat2 == j] <- mcmc_t2$batch[ss,j]
		                diag( Q_mcmc_t2_tmp ) <- 0.0 ## Make sure the diag is 0!
		                diag( Q_mcmc_t2_tmp ) <- rowSums( Q_mcmc_t2_tmp ) * -1
		                rownames(Q_mcmc_t2_tmp) <- colnames(Q_mcmc_t2_tmp) <- 1:k_t2

		                ## Another round of sampling!
		                init_nshifts <- 2000
		                n_try_t2_simmap <- 0
		                while( !"phylo" %in% class(T2[[solve_id]]) ){
		                    T2[[solve_id]] <- ratematrix::fastSimmap(tree = t2[[i]], x = y
		                                                             , pi = root_p2
		                                                             , Q = Q_mcmc_t2_tmp
		                                                             , silence = TRUE
		                                                             , max_nshifts = init_nshifts)
		                    init_nshifts <- init_nshifts * 10
		                    n_try_t2_simmap <- n_try_t2_simmap + 1
		                    if( n_try_t2_simmap > 5 ) break ## Give up!
		                }
		            }
		        }
		    }

			## Save the T1 samples.
			saveRDS(T2, file = paste0(dir_files, "/simmap_t2_rep", tag, i, ".rds") )
		}

		########################################
		## Block to compute the D quantities.
		########################################
		## Here we are using the 'obs_nrep' because of the 'read_only' parameter.

		## Need to add check for errors in the simulation replicates.
		t1_check_class <- all( sapply(T1, function(x) "phylo" %in% class(x)) )
		t2_check_class <- all( sapply(T2, function(x) "phylo" %in% class(x)) )
		if( !t1_check_class ) print( paste0("T1 trees rep ", i, " not phylo. Skipping...") )
		if( !t2_check_class ) print( paste0("T2 trees rep ", i, " not phylo. Skipping...") )

		tmp <- calcD(T1,T2)
		PD <- PD+(tmp$E_DX>=E_DX)/obs_nrep
		Pdij <- Pdij+(tmp$E_Dij>=E_Dij)/obs_nrep
		cat(".")
		## Plot at each 10 interval. Expect i to follow nrep and not obs_nrep.
		if(i %% 10 == 0 && i != nrep) cat("\n")
		dev.flush()
	}
	cat("|\nDone.\n")
	obj <- list("E(D|X)"=E_DX,"P(D)"=PD,"E(Dij)"=E_Dij,
		"P(Dij)"=Pdij)
	class(obj) <- "Dtest"
	obj
}
