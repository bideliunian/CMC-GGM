
neighbor.vtv.glasso = function(data, p, lambda){
  n = dim(data)[1]
  M = dim(data)[2]/p

  # Run ADMM and gain TPR/FPR
  G.list <- list()
  TPR.and <- rep(0,L)
  FPR.and <- rep(0,L)
  TPR.or <- rep(0,L)
  FPR.or <- rep(0,L)
  
  # default warm start ABU
  A.def <- list()
  for(j in 1:(p-1)) A.def[[j]] <- matrix(0, M, M)
  B.def <- matrix(0.1, n, M)
  U.def <- matrix(0.01, n, M)
  V <- matrix(NA, p, p)
    
  for (j in 1:p) {
    jth.range <- (((j-1)*M+1) : (j*M))
    
    Y <- data[, jth.range]
    X <- data[, -jth.range]
    A <- A.def; B <- B.def; U <- U.def
    
    grp.lasso.result <- ADMM.grplasso(X, Y, lambda,A.init=A, B.init=B, U.init=U)
    A <- grp.lasso.result$A
    B <- grp.lasso.result$B
    U <- grp.lasso.result$U
    
    # Process the estimated A into a neighbor selection vector
    A.mean <- rep(0, p-1)
    A.nonzero <- rep(0, p-1)
    for(juliet in 1:(p-1)){
      A.mean[juliet] <- mean(A[[juliet]])
      A.nonzero[juliet] <- (A.mean[juliet]!=0)
    }
    if(sum(A.nonzero)==0){ # all A are 0
      threshold <- 0
    }
    else{
      threshold <- mean(abs(A.mean[which(A.nonzero!=0)])) * 0.05
    }
    
    # Neighbor recognization
    V.j <- rep(0, p)
    for(juliet in 1:(p-1)){
      if(juliet<j){
        if(abs(A.mean[juliet]) > threshold) V.j[juliet] <- 1
      }
      else{
        if(abs(A.mean[juliet]) > threshold) V.j[juliet + 1] <- 1
      }
    }
  
    V[j,] <- V.j
  }
  
  return(V)
}
