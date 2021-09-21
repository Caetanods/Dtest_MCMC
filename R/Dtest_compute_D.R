#' Function to read the posterior predictive samples and compute the D test values.
#'
#' Note that Both T1 and T2 replicates need to be in the same directory. Note that the final test result will have sample size equal to the number of replicates in the directory. This function will not perform simulations.
#'
#' @param t1 a list of stochastic maps.
#' @param t2 a list of stochastic maps.
#' @param tag a name (pattern) used to name the files produced by the simulation.
#' @param dir_files the directory where the files will be save. You WILL need to know where this location is!
#' @param only_reps a numeric vector with the indexes of the replicates to read from the directory. If NULL, then the function will expect 'length(t1)' replicates (length(t1) == length(t2) ). Use this parameter to control which replicates will be used to compute the D test statistics.
#'
#' @return a list with the D test results.
#' @importFrom phytools getStates
#' @importFrom progress progress_bar
#' @export
Dtest_compute_D <- function(t1, t2, tag = "test", dir_files = "./partials", only_reps = NULL){

    ## This function will not make any simulation. It will only read from the directory.
    ## Function is assuming that every stochastic map replicate in the directory is good to go. We can add a step to verify the condition, but be aware that this would really increase the time to run.
    cat( "This function will read the replicates from the directory. This will take some time.\n" )

	## Perform some data checks.
	if( length(t1) != length(t2) ) stop("Length of t1 and t2 needs to be the same.")
	if( !inherits(t1,"multiPhylo") ){
	    class(t1) <- c(class(t1), "multiPhylo" )
	}
	if( !inherits(t2,"multiPhylo") ){
	    class(t2) <- c(class(t2), "multiPhylo" )
	}
	if( !inherits(t1[[1]],"simmap") ) stop("t1 needs to be a list of simmaps.")
	if( !inherits(t2[[1]],"simmap") ) stop("t2 needs to be a list of simmaps.")
	check_simmap_mat_t1 <- sapply(t1, function(x) is.matrix(x$Q) )
	check_simmap_mat_t2 <- sapply(t2, function(x) is.matrix(x$Q) )
	if( !all(check_simmap_mat_t1) ) stop("All simmaps need to have the Q matrix as the element $Q of the list. simmap format needs to be the same as in phytools::make.simmap or ratematrix::fastSimmap")
	if( !all(check_simmap_mat_t2) ) stop("All simmaps need to have the Q matrix as the element $Q of the list. simmap format needs to be the same as in phytools::make.simmap or ratematrix::fastSimmap")
	if( is.null(dir_files) ) stop('Need to provide the path to the directory with the simulation replicates.')
	if( is.null(only_reps) ){
	    run_reps <- 1:length(t1)
	} else{
	    if( !is.numeric(only_reps) ) stop("only_reps needs to NULL or a numeric vector.")
	    if( any(is.na(only_reps)) ) stop("only_reps cannot have NA values. Please check!")
	    run_reps <- only_reps
	}

	## Make a list of the files already existent in the directory.
	file_list <- list.files(path = dir_files, pattern = "*.rds", full.names = FALSE)

	## Make a list of the files for each of T1 and T2. Here we are not checking for the MCMC files.
	has_for_t1 <- vector(mode = "character", length = 0)
	has_for_t2 <- vector(mode = "character", length = 0)
	for(i in run_reps){
	    if( paste0("simmap_t1_rep", tag, i, ".rds") %in% file_list ){
	        has_for_t1 <- append(has_for_t1, paste0("simmap_t1_rep", tag, i, ".rds"))
	    }
	    if( paste0("simmap_t2_rep", tag, i, ".rds") %in% file_list ){
	        has_for_t2 <- append(has_for_t2, paste0("simmap_t2_rep", tag, i, ".rds"))
	    }
	}
	## Now we will use the minimum sample size. Note that we need to match for T1 and for T2, but we don't need match between T1 and T2.
	global_ss <- min( length(has_for_t1), length(has_for_t2) )
	## Now we need to create a vector of file names with the size of "global_ss". This is going to be the number of replicates used here.
	has_for_t1 <- has_for_t1[1:global_ss]
	has_for_t2 <- has_for_t2[1:global_ss]

	## Get the number of states on each of the vectors.
	levs1 <- sort(unique(as.vector(getStates(t1,"tips"))))
	k1 <- length(levs1)
	levs2 <- sort(unique(as.vector(getStates(t2,"tips"))))
	k2 <- length(levs2)

	## Calculates the D metric given the observed stochastic maps.
	obj <- calcD(t1,t2) ## Note that we have a distribution of stochastic maps here.
	E_DX <- obj$E_DX
	E_Dij <- obj$E_Dij

	## posterior prediction (Here is where we do the simulations)
	PD <- 0
	## Matrix to store the results.
	Pdij <- matrix(0,k1,k2,dimnames=list(levs1,levs2))

	## The number of replicates is the one we can read from the directory:
	obs_nrep <- global_ss

	## Create the progress bar:
	pb <- progress_bar$new(total = global_ss)

	for( i in 1:global_ss){
	    ## Both T1 and T2 replicates need to be in the same directory.
	    T1 <- readRDS( file.path(dir_files, has_for_t1[i]) )
	    T2 <- readRDS( file.path(dir_files, has_for_t2[i]) )
	    tmp <- calcD(T1,T2)
	    PD <- PD+(tmp$E_DX>=E_DX)/obs_nrep
	    Pdij <- Pdij+(tmp$E_Dij>=E_Dij)/obs_nrep
	    ## Update the progress bar:
	    pb$tick()
	}

	obj <- list("E(D|X)"=E_DX, "P(D)"=PD, "E(Dij)"=E_Dij, "P(Dij)"=Pdij, sample_size=obs_nrep)
	class(obj) <- "Dtest_stats"
	return(obj)
}
