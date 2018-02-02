#' Select the most representative transmission tree from a MCMC output with multiple sperate trees after running inferMultiTTree
#' @param record Output from inferMultiTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return The indexes of the selected transmission trees
#' @export 
selectMultiTTree = function(record,burnin=0.5)
{
  n.trees <- length(record[[1]]$ctree.list)
  
  indexes <- lapply(1:n.trees, function(i){
    for (j in 1:length(record)){
      record[[j]]$ctree <- record[[j]]$ctree.list[[i]]
    }
    TransPhylo::selectTTree(record,burnin=burnin)
  })
  names(indexes) <- 1:n.trees
  return(indexes)
}

  