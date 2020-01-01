# generates divisions of the data for cross validation
# n: number of observations
# K: number of folds
# nrep: how often should the CV procedure be repeated
# value:
# list with one entry for every nrep
# every nrep entry is a list with K entries
# every K entry is a vector of the length n with 1 for the observations in the
# training set and 0 for the observations in the test set
makeCVdivision <- function(n, K = 5, nrep = 3) {
    
    ngroup <- rep(floor(n/K), times=K)
    if(n - sum(ngroup) != 0)
      ngroup[1:(n - sum(ngroup))] <- ngroup[1:(n - sum(ngroup))] + 1
    
    cvlist <- list()
    
    for(i in 1:nrep) {
      indgroup <- sample(rep(1:K, times=ngroup))
      cvlist[[i]] <- list()
      for(j in 1:K) {
        tempvec <- rep(0, n)
        tempvec[indgroup!=j] <- 1   
        cvlist[[i]][[j]] <- tempvec
      }
    }
  
  return(cvlist)

}
