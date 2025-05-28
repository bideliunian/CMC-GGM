graph.est <- function(S, n, p, sparsity = NULL){
  ## score: fpca score vectors with dimension n by (p*M)
  ## method: {'glasso','thresholding'}
  ## trans: {'cmcopula', 'copula', 'gaussian'}
  ## L: number of lambdas
  M = nrow(S) / p
  
  tol.Theta <- 1e-2
  tol.beta <- 1e-3
  maxit.Theta <- 10
  maxit.beta <- 10
  tol.beta_k <- 1e-7
  maxit.beta_k <- 1500
  method <- 1
  L = 8
  
  lambda = seq(0.2, 0.5, length = L)
  
  rho = eigen(S)$values[1]
  
  cl <- makeCluster(L)
  registerDoParallel(cl)
  
  result <- foreach(l = 1:L, .combine="c", .packages=c("matrixcalc","igraph","BB")) %dopar% {
    func.path <- "D:/Research/VcopularGGM/Codes/Functions"
    source(paste(func.path,"ns_lasso.R", sep="/"))
    bglasso.list <- bglasso(S = S, lambda = lambda[l], n = n, p = p,
                            tol.Theta=tol.Theta, tol.beta=tol.beta,
                            tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                            maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)
    list(bglasso.list)
  }
  
  stopCluster(cl)
  
  bic.list <- lapply(result, "[[", 'bic')
  aic.list <- lapply(result, "[[", 'aic')
  est.support.list <- lapply(result, "[[", 'alpha')
  Theta.list <- lapply(result, "[[", 'Theta')
  iter.list <- lapply(result, "[[", 'iter.Theta')
  
  print(paste("aic:", as.vector(unlist(aic.list))))
  print(paste("bic:", as.vector(unlist(bic.list))))
  print(paste("iter:", as.vector(unlist(iter.list))))
  
  if(!is.null(sparsity)){
    sparsity.list <- lapply(est.support.list, function(x) sum(x[upper.tri(x)] > 0)/(p*(p-1)/2))
    sparsity.list <- as.vector(unlist(sparsity.list))
    print(sparsity.list)
    index <- which.min(abs(sparsity.list - sparsity))
    sparsity.graph <- sparsity.list[index]
    est.support <- est.support.list[[index]] 
  }
  else{
    bic.list <- as.vector(unlist(bic.list))
    print(bic.list)
    index <- which.min(bic.list)
    est.support <- est.support.list[[index]]
    sparsity.graph <- sum(est.support[upper.tri(est.support)] > 0)/(p*(p-1)/2)
  }
  
  est.Theta <- Theta.list[[index]] 
  
  list(est.support=est.support, sparsity=sparsity.graph, Theta = est.Theta)
}