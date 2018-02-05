#' Build a matrix of probability of who infected whom from a MCMC output
#' @param record Output from inferMultiTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return Matrix of probability of who infected whom
#' @export
computeMatWIWMulti = function(record,burnin=0.5)
{
  n.trees <- length(record[[1]]$ctree.list)
  
  mat.wiws <- lapply(1:n.trees, function(i){
    for (j in 1:length(record)){
      record[[j]]$ctree <- record[[j]]$ctree.list[[i]]
      record[[j]]$source <- record[[j]]$source.list[[i]]
    }
    TransPhylo::computeMatWIW(record, burnin=burnin)
  })
  names(mat.wiws) <- 1:n.trees
  
  isolates <- unlist(lapply(mat.wiws, rownames))
  mat.wit <- matrix(Inf, nrow = length(isolates), ncol = length(isolates), 
                    dimnames = list(rownames=isolates, colnames=isolates))
  
  for (i in 1:n.trees){
    mat.wit[cbind(rownames(mat.wiws[[i]])[row(mat.wiws[[i]])], 
                  colnames(mat.wiws[[i]])[col(mat.wiws[[i]])])] <- c(mat.wiws[[i]])
  }
  return(mat.wit)
}

#' Build a matrix indicating for each pairs of individuals how many intermediates there are in the transmission chain
#' @param record Output from inferMultiTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return Matrix of intermediates in transmission chains between pairs of hosts
#' @export
computeMatTDistMulti = function(record,burnin=0.5)
{
  n.trees <- length(record[[1]]$ctree.list)
  
  mat.tdists <- lapply(1:n.trees, function(i){
    for (j in 1:length(record)){
      record[[j]]$ctree <- record[[j]]$ctree.list[[i]]
      record[[j]]$source <- record[[j]]$source.list[[i]]
    }
    TransPhylo::computeMatWIW(record, burnin=burnin)
  })
  names(mat.tdists) <- 1:n.trees
  
  isolates <- unlist(lapply(mat.tdists, rownames))
  mat.tdist <- matrix(Inf, nrow = length(isolates), ncol = length(isolates), 
                    dimnames = list(rownames=isolates, colnames=isolates))
  
  for (i in 1:n.trees){
    mat.tdist[cbind(rownames(mat.tdists[[i]])[row(mat.tdists[[i]])], 
                  colnames(mat.tdists[[i]])[col(mat.tdists[[i]])])] <- c(mat.tdists[[i]])
  }
  return(mat.tdist)
}

