#' @importFrom phytools sim.Mk
reject_sample_MK <- function(tree, Q, root_pi){
    ## Performs MK simulation and returns only if data has the same number of states as the generating matrix.
    n_target <- ncol(Q)
    another_sim <- TRUE
    n_trials <- 0
    while( another_sim ){
        ss <- sim.Mk(tree = tree, Q = Q, anc = sample(x = colnames(Q), size = 1, prob = root_pi))
        another_sim <- length(unique(ss)) != n_target
        n_trials <- n_trials + 1
        ## print( paste0("Ntrial:", n_trials, " - Ntraits: ", length(unique(ss))) )
        if( n_trials > 1000 ){
            another_sim <- FALSE
            ## warning( "Unable to simulate n traits!" )
            print( "Unable to simulate n traits!" )
        }
    }
    return( ss )
}
