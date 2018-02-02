#' Infer transmission tree given multiple phylogenetic trees
#' @param ptree.list List of phylogenetic tree
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
#' @param update.w.scale Whether of not to update the parameter w.scale
#' @param update.w.shape Whether of not to update the parameter w.shape
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param startNeg Starting value of within-host coalescent parameter Ne*g
#' @param startOff.r Starting value of parameter off.r
#' @param startOff.p Starting value of parameter off.p
#' @param startPi Starting value of sampling proportion pi
#' @param updateNeg Whether of not to update the parameter Ne*g
#' @param updateOff.r Whether or not to update the parameter off.r
#' @param updateOff.p Whether or not to update the parameter off.p
#' @param updatePi Whether or not to update the parameter pi
#' @param updateTTree Whether or not to update the transmission tree
#' @param optiStart Whether or not to optimise the MCMC start point
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @param delta.t Interval to calculate integral.
#' @return posterior sample set of transmission trees
#' @export

inferMultiTTree = function(ptree.list, w.shape=2, w.scale=1, ws.shape=w.shape, ws.scale=w.scale,
                      mcmcIterations=1000, thinning=1, startNeg=100/365, startOff.r=1,
                      startOff.p=0.5, startPi=0.5, updateNeg=TRUE, updateOff.r=TRUE,
                      updateOff.p=FALSE, updatePi=TRUE, updateTTree=TRUE,
                      update.w.scale=FALSE, update.w.shape=FALSE, optiStart=TRUE, dateT=Inf,
                      delta.t=0.01) {
  
  if (update.w.scale && updatePi) stop("Error: can not estimate pi and w.scale simultaneously.")
  
  for (ptree in ptree.list){
    ptree$ptree[,1]=ptree$ptree[,1]+runif(nrow(ptree$ptree))*1e-10#Ensure that all leaves have unique times
    for (i in (ceiling(nrow(ptree$ptree)/2)+1):nrow(ptree$ptree)) for (j in 2:3) 
      if (ptree$ptree[ptree$ptree[i,j],1]-ptree$ptree[i,1]<0) 
        stop("The phylogenetic tree contains negative branch lengths!")
  }
  
  #MCMC algorithm
  neg <- startNeg
  off.r <- startOff.r
  off.p <- startOff.p
  pi <- startPi
  ctree.list <- lapply(ptree.list, function(ptree){
    makeCtreeFromPTree(ptree,ifelse(optiStart,off.r,NA),off.p,neg,pi,w.shape,w.scale,ws.shape,ws.scale,dateT)#Starting point 
  })

  ttree.list <- lapply(ctree.list, extractTTree)
  record <- vector('list',mcmcIterations/thinning)
  pTTree.list <- lapply(ttree.list, function(ttree) {
    probTTree(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT,delta.t)
  })
  pPTree.list <- lapply(ctree.list, probPTreeGivenTTreectree,neg)
  pb <- txtProgressBar(min=0,max=mcmcIterations,style = 3)
  for (i in 1:mcmcIterations) {#Main MCMC loop
    if (i%%thinning == 0) {
      #Record things 
      setTxtProgressBar(pb, i)
      #message(sprintf('it=%d,neg=%f,off.r=%f,off.p=%f,pi=%f,Prior=%e,Likelihood=%f,n=%d',i,neg,off.r,off.p,pi,pTTree,pPTree,nrow(extractTTree(ctree))))
      record[[i/thinning]]$ctree.list <- ctree.list
      record[[i/thinning]]$pTTree <- sum(pTTree.list)
      record[[i/thinning]]$pPTree <- sum(pPTree.list)
      record[[i/thinning]]$neg <- neg 
      record[[i/thinning]]$off.r <- off.r
      record[[i/thinning]]$off.p <- off.p
      record[[i/thinning]]$pi <- pi
      record[[i/thinning]]$w.shape <- w.shape
      record[[i/thinning]]$w.scale <- w.scale
      record[[i/thinning]]$ws.shape <- ws.shape
      record[[i/thinning]]$ws.scale <- ws.scale
      record[[i/thinning]]$source <- lapply(ctree.list, function(ctree) {
        source <- ctree$ctree[ctree$ctree[which(ctree$ctree[,4]==0),2],4]
        if (source<=length(ctree$nam)) {
          source=ctree$nam[source] 
        } else { 
          source='Unsampled'
        }
        return(source)
      })
    }
    
    if (updateTTree) {
      tree.to.move <- sample.int(length(ctree.list),1)
      #Metropolis update for transmission tree
      prop <- proposal(ctree.list[[tree.to.move]]$ctree) 
      ctree2 <- list(ctree=prop$tree,nam=ctree.list[[tree.to.move]]$nam)
      ttree2 <- extractTTree(ctree2)
      pTTree2 <- probTTree(ttree2$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT,delta.t) 
      pPTree2 <- probPTreeGivenTTree(ctree2,neg) 
      if (log(runif(1)) < log(prop$qr)+ pTTree2 + pPTree2 - pTTree.list[tree.to.move] - pPTree.list[tree.to.move])  { 
        ctree.list[tree.to.move] <- ctree2
        ttree.list[tree.to.move] <- ttree2
        pTTree.list[tree.to.move] <- pTTree2 
        pPTree.list[tree.to.move] <- pPTree2 
      } 
    }
    
    if (updateNeg) {
      #Metropolis update for Ne*g, assuming Exp(1) prior 
      neg2 <- abs(neg + (runif(1)-0.5)*0.5)
      pPTree2.list <- lapply(ctree.list, probPTreeGivenTTree, neg2)
      if (log(runif(1)) < sum(pPTree2.list)-sum(pPTree.list)-neg2+neg)  {
        neg <- neg2
        pPTree.list <- pPTree2.list
        } 
    }
    
    if (update.w.scale) {
      #Metropolis update for w.scale, assuming Exp(1) prior shifted by 1e-3 or just under a day 
      w.scale.2 <- abs(w.scale + (runif(1)-0.5)*0.5)
      if (w.scale.2<0.01) w.scale.2=0.02-w.scale.2
      pTTree2.list <- lapply(ttree.list, function(ttree) {
        probTTree(ttree$ttree,off.r,off.p,pi,w.shape,w.scale.2,ws.shape,ws.scale,dateT,delta.t)
      })
      if (log(runif(1)) < sum(pTTree2.list)-sum(pTTree.list)+shifted_gamma_prior(w.scale.2)-shifted_gamma_prior(w.scale)) {
        w.scale <- w.scale.2
        pTTree.list <- pTTree2.list
      }
    }
    
    if (update.w.shape) {
      #Metropolis update for w.shape, assuming Exp(1) prior shifted by 1e-3 or just under a day 
      w.shape.2 <- abs(w.shape + (runif(1)-0.5)*0.5)
      if (w.shape.2<0.01) w.shape.2=0.02-w.shape.2
      pTTree2.list <- lapply(ttree.list, function(ttree){
        probTTree(ttree$ttree,off.r,off.p,pi,w.shape.2,w.scale,ws.shape,ws.scale,dateT,delta.t)
      })
      if (log(runif(1)) < sum(pTTree2.list)-sum(pTTree.list)+shifted_gamma_prior(w.shape.2)-shifted_gamma_prior(w.shape)) {
        w.shape <- w.shape.2
        pTTree.list <- pTTree2.list
      }
    }
    
    if (updateOff.r) {
      #Metropolis update for off.r, assuming prior exp(1)  
      off.r2 <- abs(off.r + (runif(1)-0.5)*0.5)
      pTTree2.list <- lapply(ttree.list, function(ttree){
        probTTree(ttree$ttree,off.r2,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT,delta.t)
      })
      if (log(runif(1)) < sum(pTTree2.list)-sum(pTTree.list)-off.r2+off.r)  {
        off.r <- off.r2
        pTTree.list <- pTTree2.list
        }
    }
    
    if (updateOff.p) {
      #Metropolis update for off.p, assuming Unif(0,1) prior 
      off.p2 <- abs(off.p + (runif(1)-0.5)*0.1)
      if (off.p2>1) off.p2=2-off.p2
      pTTree2.list <- lapply(ttree.list, function(ttree) {
        probTTree(ttree$ttree,off.r,off.p2,pi,w.shape,w.scale,ws.shape,ws.scale,dateT,delta.t)
      })
      if (log(runif(1)) < pTTree2-pTTree)  {
        off.p <- off.p2
        pTTree.list <- pTTree2.list
        }
    }
    
    if (updatePi) {
      #Metropolis update for pi, assuming Unif(0.01,1) prior 
      pi2 <- pi + (runif(1)-0.5)*0.1
      if (pi2<0.01) pi2=0.02-pi2
      if (pi2>1) pi2=2-pi2
      pTTree2.list <- lapply(ttree.list, function(ttree){
        probTTree(ttree$ttree,off.r,off.p,pi2,w.shape,w.scale,ws.shape,ws.scale,dateT,delta.t)
      })
      if (log(runif(1)) < sum(pTTree2.list)-sum(pTTree.list))  {
        pi <- pi2
        pTTree.list <- pTTree2.list
        }
    }
    
  }#End of main MCMC loop
  
  return(record)
}

shifted_exp_prior <- function(x){
  L <- 1e-3
  return(-(x-L))
}

shifted_gamma_prior <- function(x){
  L <- 1e-2
  shape <- 1
  scale <- 2
  return((shape-1)*log(x-L)-(x-L)/scale-shape*log(scale)-log(gamma(shape)))
}
