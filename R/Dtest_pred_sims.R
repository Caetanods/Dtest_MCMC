#' Function to perform the MCMC posterior predictive simulations
#'
#' This function will only run the simulations. After running the posterior predictive simulations, another function will be used to compute the test statistics. Note that the function will use rejection sampling to make sure that each predictive simulation has the same number of states as the original data (maybe we can relax this assumption and add a parameter to control this).
#'
#' @param t1 a list of stochastic maps. Length of t1 and t2 need to be the same.
#' @param t2 a list of stochastic maps. Length of t1 and t2 need to be the same.
#' @param tag a name (pattern) used to name the files produced by the simulation.
#' @param dir_files the directory where the files will be save. You WILL need to know where this location is!
#' @param only_reps a numeric vector with the indexes of the replicates to run. If NULL, then the function will do a number of replicates equal to the length of t1 ( length(t1) == length(t2) ). Use this parameter to spread replicates through multiple processors or computers.
#' @param ngen number of generations for each MCMC chain.
#' @param nsim number of stochastic mapping samples from each replicate.
#' @param model_mat1 the model for the trait 1. Same format as in the "model" parameter for the function "ace" in the package "ace".
#' @param model_mat2 the model for the trait 2. Same format as in the "model" parameter for the function "ace" in the package "ace".
#' @param root_p1 the vector of root probabilities for trait 1.
#' @param root_p2 the vector of root probabilities for trait 2.
#' @param prior_exp the parameter for the exponential prior for the rates of transition of the Markov model fitted to the data.
#' @param scale_mcmc the scaling factor for the Metropolis step of the MCMC. See function "metrop" in package "mcmc".
#'
#' @return Function will write a large number of ".rds" files to the directory. Please make sure that you have read and write permissions to the "dir_files" location! Also make sure to record the location of the "dir_files" directory. It also returns a list with information on replicates with problems.
#' @importFrom ratematrix mergeSimmap
#' @importFrom ratematrix fastSimmap
#' @importFrom corHMM corHMM
#' @importFrom mcmc metrop
#' @importFrom phytools getStates
#' @importFrom progress progress_bar
#' @export
Dtest_pred_sims <- function(t1, t2, tag = "test", dir_files, ngen = 5000, nsim = 100, only_reps = NULL, model_mat1, model_mat2, root_p1, root_p2, prior_exp = 0.5, scale_mcmc = 1){

    cat( "This function will make posterior predictive simulations. This will take some time.\n" )

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
    if( is.null(only_reps) ){
        run_reps <- 1:length(t1)
    } else{
        if( !is.numeric(only_reps) ) stop("only_reps needs to NULL or a numeric vector.")
        if( any(is.na(only_reps)) ) stop("only_reps cannot have NA values. Please check!")
        run_reps <- only_reps
    }
    if( is.null(dir_files) ) stop('Need to provide a path to write MCMC samples. Use dir_files="." to write files to current directory.')
    ## Need to extend it to accept the "madftiz" option.
    if( is.numeric( root_p1 ) ){
        if( length(root_p1) != ncol(model_mat1) ) stop("Length of root_p1 needs to be equal to dimension of model_mat1")
    }
    if( is.numeric( root_p2 ) ){
        if( length(root_p2) != ncol(model_mat2) ) stop("Length of root_p1 needs to be equal to dimension of model_mat1")
    }

    ## Creates a directory to store the partial results.
    dir.create( dir_files, showWarnings = FALSE )

    ## Get the number of states on each of the vectors.
    levs1 <- sort(unique(as.vector(getStates(t1,"tips"))))
    k1 <- length(levs1)
    levs2 <- sort(unique(as.vector(getStates(t2,"tips"))))
    k2 <- length(levs2)

    ## Get the number of parameters for each of the models:
    npar_t1 <- sum(unique(c(model_mat1)) > 0, na.rm = TRUE)
    npar_t2 <- sum(unique(c(model_mat2)) > 0, na.rm = TRUE)

    ## Create a list to store information about the replicates that did not work.
    t1_error_list <- list()
    t2_error_list <- list()

    ## Create the progress bar:
    pb <- progress_bar$new(total = length(run_reps))

    ## For each of the posterior predictive reps we need to do a full MCMC.
    for(i in run_reps){

        ## Need to fix the rownames of the Q matrix:
        Q_t1 <- t1[[i]]$Q
        k_t1 <- ncol( Q_t1 )
        rownames( Q_t1 ) <- colnames( Q_t1 )
        Q_t2 <- t2[[i]]$Q
        k_t2 <- ncol( Q_t2 )
        rownames( Q_t2 ) <- colnames( Q_t2 )

        ## Simulation under the model to make data for a posterior predictive. (For t1 and t2)

        ## Instead of making a simple simulation, here we are using rejection sampling.
        ## y <- phytools::sim.Mk(t2[[i]], Q_t2, anc = sample(x = colnames(Q_t2), size = 1, prob = root_p2))
        x <- reject_sample_MK(tree = t1[[i]], Q = Q_t1, root_pi = root_p1)
        x <- setNames(as.character(x), names(x))
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

        ## Run the MCMC and save the rds file.
        ## Get the MLE as starting state.
        t1_tmp <- mergeSimmap(phy = t1[[i]], drop.regimes = TRUE)
        invisible(capture.output( start_mle_guess_corHMM <- corHMM(phy = t1_tmp, data = data.frame(spp = names(x_num), trait = as.numeric(x_num)), rate.cat = 1, rate.mat = model_mat1, root.p = root_p1) ))
        start_mle_guess <- sapply(1:max(model_mat1, na.rm = TRUE), function(x) start_mle_guess_corHMM$solution[which(model_mat1 == x)][1])
        if( any(start_mle_guess < 1e-6) ){
            ## Bounce the starting state to a resonable value if it is zero.
            start_mle_guess[start_mle_guess < 1e-6] <- 1e-6
        }

        ## Run the MCMC.
        mcmc_t1 <- metrop(obj = post_fn_t1, initial = log(start_mle_guess), nbatch = ngen
                          , blen = 1, nspac = 1, scale = scale_mcmc)
        mcmc_t1$batch <- exp(mcmc_t1$batch) ## Need to untransform the parameter space.
        ## Save the partial result.
        ## Need to remember that this posterior is in a distinct format.
        saveRDS(mcmc_t1, file = paste0(dir_files, "/mcmc_t1_rep", tag, i, ".rds") )

        ## Sample Q matrices from the MCMC in order to make the stochastic maps.
        ss <- sample(x = seq(from = round(ngen/2), to = ngen, by = 1), size = nsim, replace = FALSE)
        dt_t1 <- mcmc_t1$batch[ss,]
        ## If model is ER, then dt_t1 will be a vector and not a matrix:
        if( is.vector(dt_t1) && is.atomic(dt_t1) ){
            ## Transform the sample in a list of Q matrices.
            Q_mcmc_t1 <- list()
            for( w in 1:nsim ){
                ## Need to reconstruct using the model matrix.
                mm <- model_mat1
                ## ER model only!
                mm[which(model_mat1 == 1)] <- dt_t1[i]
                diag( mm ) <- 0.0 ## Make sure the diag is 0!
                diag( mm ) <- rowSums( mm ) * -1
                rownames(mm) <- colnames(mm) <- 1:k_t1
                Q_mcmc_t1[[w]] <- mm
            }
        } else{
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
        }

        ## Define a recursive function to run the stochastic maps. If fastSimmap runs into the limit, then ti will return NA.
        run_fast_simmap_t1 <- function(q_mat, limit_shifts = 200){
            ## If the limit becames too high, give up.
            if( limit_shifts > 20000 ) return( NA )
            ss_sim <- fastSimmap(tree = t1[[i]], x = x, pi = root_p1, Q = q_mat
                                 , silence = TRUE, max_nshifts = limit_shifts)
            if( "phylo" %in% class(ss_sim) ){
                return( ss_sim )
            } else{
                ## Try again:
                run_fast_simmap_t1(q_mat = q_mat, limit_shifts = limit_shifts * 10)
            }
        }

        T1 <- lapply(Q_mcmc_t1, function(xxx) run_fast_simmap_t1(q_mat = xxx))
        ## Check for problems in the stochastic map simulation:
        t1_check_class <- sapply(T1, function(x) "phylo" %in% class(x))
        if( !all(t1_check_class) ){
            failed_t1 <- which(!t1_check_class)
            cat("For replicate", i, "these t1 samples failed:", failed_t1, sep = " ")
            print( "Transition matrix for t1 might have rates that are too fast. Check MCMC samples." )
        }

        ## Save only if all replicates worked. Later we can introduce a way to deal with the NAs.
        if( all(t1_check_class) ){
            t1_error_list <- append(t1_error_list, NA)
            saveRDS(T1, file = paste0(dir_files, "/simmap_t1_rep", tag, i, ".rds") )
        } else{
            t1_error_list <- append(t1_error_list, which(!t1_check_class))
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

        ## Run the MCMC and save the rds file.
        ## Get the MLE as starting state.
        t2_tmp <-  mergeSimmap(phy = t2[[i]], drop.regimes = TRUE)
        invisible(capture.output( start_mle_guess_corHMM <- corHMM(phy = t2_tmp, data = data.frame(spp = names(y_num), trait = as.numeric(y_num)), rate.cat = 1, rate.mat = model_mat2, root.p = root_p2) ))
        start_mle_guess <- sapply(1:max(model_mat2, na.rm = TRUE), function(x) start_mle_guess_corHMM$solution[which(model_mat2 == x)][1])
        if( any(start_mle_guess < 1e-6) ){
            start_mle_guess[start_mle_guess < 1e-6] <- 1e-6
        }

        ## Run the MCMC.
        mcmc_t2 <- metrop(obj = post_fn_t2, initial = log(start_mle_guess), nbatch = ngen
                          , blen = 1, nspac = 1, scale = scale_mcmc)
        mcmc_t2$batch <- exp(mcmc_t2$batch) ## Need to unstrasform the MCMC samples
        ## Save the partial result.
        saveRDS(mcmc_t2, file = paste0(dir_files, "/mcmc_t2_rep", tag, i, ".rds") )

        ss <- sample(x = seq(from = round(ngen/2), to = ngen, by = 1), size = nsim, replace = FALSE)
        dt_t2 <- mcmc_t2$batch[ss,]
        ## If model is ER, then dt_t2 will be a vector and not a matrix:
        if( is.vector(dt_t2) && is.atomic(dt_t2) ){
            ## Transform the sample in a list of Q matrices.
            Q_mcmc_t2 <- list()
            for( w in 1:nsim ){
                ## Need to reconstruct using the model matrix.
                mm <- model_mat2
                ## ER model only!
                mm[which(model_mat2 == 1)] <- dt_t2[i]
                diag( mm ) <- 0.0 ## Make sure the diag is 0!
                diag( mm ) <- rowSums( mm ) * -1
                rownames(mm) <- colnames(mm) <- 1:k_t2
                Q_mcmc_t2[[w]] <- mm
            }
        } else{
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
        }

        ## Make the stochastic maps and save the rds file.
        run_fast_simmap_t2 <- function(q_mat, limit_shifts = 200){
            ## If the limit becames too high, give up.
            if( limit_shifts > 20000 ) return( NA )
            ss_sim <- fastSimmap(tree = t2[[i]], x = x, pi = root_p2, Q = q_mat
                                 , silence = TRUE, max_nshifts = limit_shifts)
            if( "phylo" %in% class(ss_sim) ){
                return( ss_sim )
            } else{
                ## Try again:
                run_fast_simmap_t2(q_mat = q_mat, limit_shifts = limit_shifts * 10)
            }
        }

        T2 <- lapply(Q_mcmc_t2, function(xxx) run_fast_simmap_t2(q_mat = xxx))
        ## Check for problems in the stochastic map simulation:
        t2_check_class <- sapply(T2, function(x) "phylo" %in% class(x))
        if( !all(t2_check_class) ){
            failed_t2 <- which(!t2_check_class)
            cat("For replicate", i, "these t2 samples failed:", failed_t2, sep = " ")
            print( "Transition matrix for t2 might have rates that are too fast. Check MCMC samples." )
        }

        ## Save the T2 samples only if they all worked. Otherwise skip.
        if( all(t2_check_class) ){
            t2_error_list <- append(t2_error_list, NA)
            saveRDS(T2, file = paste0(dir_files, "/simmap_t2_rep", tag, i, ".rds") )
        } else{
            t2_error_list <- append(t2_error_list, which(!t2_check_class))
        }

        ## Add one tick to the progress bar:
        pb$tick()
    }
    return( list( failed_t1 = t1_error_list, failed_t2 = t2_error_list) )
}
