##############################################
#   Function to calculate FPC score matrix   #
#   from observation 3D array h              #
#   based on FPCA basis using j'th node      #
##############################################
# To select neighborhood for each j, recalculate an FPC score matrix


FPC.score.j <- function(h, j, M=5){
  # h: 3D array h[i,j,k]
  # j: node index
  # M: number of basis
  
  n <- dim(h)[1]
  p <- dim(h)[2]
  tau <- dim(h)[3]
  obs.time <- seq(1/tau, 1, 1/tau)
  
  # Step 1. gain FPCA basis from h[,j,]
  obs.val.matrix.j <- matrix(0, nrow=tau, ncol=n)
  for (i in 1:n){
    obs.val.matrix.j[, i] <- as.vector(h[i, j, ])
  }
  bspline.basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=M)
  # create functional data object
  fd.object.array.j <- Data2fd(argvals=obs.time, y=obs.val.matrix.j, basisobj=bspline.basis)
  # get fpca basis for h[,j,]
  fpca.basis.j <- pca.fd(fd.object.array.j, nharm=M)$harmonics
  # this fpca basis contains M functions
  
  # Step 2. convert h to functional data object
  # Step 3. calculate inner product on the IDENTICAL fpca.basis.j (different from traditional fpca)
  X.j <- matrix(NA, nrow=n, ncol=(p-1)*M)
  for(juliet in 1:p){
    obs.val.matrix.juliet <- matrix(0, nrow=tau, ncol=n)
    for(i in 1:n){
      obs.val.matrix.juliet[, i] <- as.vector(h[i, juliet, ])
    }
    fd.object.array.juliet <- Data2fd(argvals=obs.time, y=obs.val.matrix.juliet, basisobj=bspline.basis)
    # the fd object above contains n functions
    
    inn.prod.juliet <- inprod(fd.object.array.juliet, fpca.basis.j) # a n*M matrix
    
    if(juliet == j){
      Y.j <- inn.prod.juliet
    }else if(juliet < j){
      juliet.range <- (((juliet-1)*M+1) : (juliet*M))
      X.j[ , juliet.range] <- inn.prod.juliet
    }else{
      juliet.range <- (((juliet-2)*M+1) : ((juliet-1)*M))
      X.j[ , juliet.range] <- inn.prod.juliet
    }
  }
  
  return(list(Y=Y.j, X=X.j))
}

blockwise_Frob <- function(ma, M){
  # ma: matrix
  # block size
  p <- dim(ma)[1]
  if (p %% M != 0){
    stop("dimension of matrix cannot be divided by block size")
  }
  n.b <- p / M
  result <- matrix(0, nrow=n.b, ncol=n.b)
  for (i in 1:n.b){
    for (j in 1:n.b){
      result[i, j] <- norm(ma[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)], "F")
    }
  }
  return(result)
}

object_func <- function(S, Theta, gamma, M){
  
  part1 <- -log(det(Theta))
  part2 <- sum(diag(S %*% Theta))
  
  block_frob <- blockwise_Frob(Theta, M)
  part3 <- gamma * (sum(block_frob) - sum(diag(block_frob)))
  
  return(part1+part2+part3)
}

ProxAlg_FGM <- function(S, p, M, gamma, Eta="Auto", n.iteration=2000){
  # Optimzation function to solve the optimization problem in FGM
  # Use Proximal Algorithm
  # 
  # Args:
  #   S: the estimated covariance matrix of Graph
  #   p: the number of vertices
  #   M: the number of principle components
  #   gamma: hyperparameter.
  #   Eta: hyperparamter. If Eta="Auto", then Eta would be set as inverse of the 
  #        product of the maximum singular value of cov.X and cov.Y
  #   n.iteration: maximum times of iteration
  #
  # Returns:
  #   A list:
  #     ThetaMathat: The estimated differential matrix between ThetaX and ThetaY
  #     blockFrob: A matrix of Frobenius norm of each block matrix
  #     Support: Support Matrix
  if(Eta == "Auto"){
    Eta <- 0.01 # this needs more careful design
  }
  
  v.obj <- c()
  
  converge.indicator <- FALSE
  
  Theta.old <- diag(1, nrow=(p*M), ncol=(p*M))
  for (t in 1:n.iteration){
    obj.old <- object_func(S, Theta.old, gamma, M)
    v.obj <- c(v.obj, obj.old)
    
    Theta <- Theta.old
    
    # calculate gradient
    grad.Mat <- S - round(solve(Theta), 10)  # round to solve asymmetric problem due to 
    # numerical issue
    A <- Theta - Eta * grad.Mat
    
    # update Theta
    for (i in 1:p){
      for (j in 1:p){
        Asub <- A[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]
        Fnorm <- norm(Asub,"F")
        if (i != j){
          if(Fnorm <= gamma*Eta){
            Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, nrow=M, ncol=M)
          } else {
            Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- round(((Fnorm - gamma*Eta) / Fnorm) * Asub, 9)
          }
        } else {
          Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- round(Asub, 9)
        }
      }
    }
    
    Theta.new <- Theta
    obj.new <- object_func(S, Theta.new, gamma, M)
    
    if(is.na(obj.new)| obj.new == -Inf){
      Theta.new <- Theta.old
      Theta.hat <- 0.5 * (Theta.new + t(Theta.new))
      Theta.frob <- blockwise_Frob(Theta.hat, M)
      
      SupMat <- matrix(0, nrow=p, ncol=p)
      SupMat[Theta.frob > 0] <- 1
      
      return (list(ThetaMathat=Theta.hat, blockFrob=Theta.frob, Support=SupMat, 
                   converge=converge.indicator, num.iter=t))
      break
    }
    
    if (abs((obj.new-obj.old)/(obj.old+1e-15)) < 1e-3){
      converge.indicator <- TRUE
      Theta.hat <- 0.5 * (Theta.new + t(Theta.new))
      Theta.frob <- blockwise_Frob(Theta.hat, M)
      
      SupMat <- matrix(0, nrow=p, ncol=p)
      SupMat[Theta.frob > 0] <- 1
      
      return (list(ThetaMathat=Theta.hat, blockFrob=Theta.frob, Support=SupMat, 
                   converge=converge.indicator, num.iter=t))
      break
    }
    Theta.old <- Theta.new
  }
  
  Theta.hat <- 0.5 * (Theta.new + t(Theta.new))
  Theta.frob <- blockwise_Frob(Theta.hat, M)
  
  SupMat <- matrix(0, nrow=p, ncol=p)
  SupMat[Theta.frob > 0] <- 1
  
  return (list(ThetaMathat=Theta.hat, blockFrob=Theta.frob, Support=SupMat, 
               converge=converge.indicator, num.iter=t))
}


# a function to calculate TP, FP, TN, FN
# and return precision, recall, TPR and FPR.
# Input:
#   G.true, the true p*p adjacency matrix
#   G.mat, the estimated p*p adjacency matrix
#   type,
#     AND: when two nodes both recognize each other as neighbors
#     OR: when either of them recognize each other as neighbors
# Output:
#   prec, precision
#   rec, recall
#   TPR and FPR, as their names

prec.rec <- function(G.true, G.mat, type=c("AND","OR")){
  
  p <- nrow(G.true)
  TP <- 0; TN <- 0; FP <- 0; FN <- 0
  
  if(type=="AND"){
    for(i in 1:p){
      for(j in 1:p){
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==1 & G.mat[j,i]==1)) TP <- TP + 1
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==0 & G.mat[j,i]==0)) FN <- FN + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==1 & G.mat[j,i]==1)) FP <- FP + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==0 & G.mat[j,i]==0)) TN <- TN + 1
      }
    }
  }else if(type=="OR"){
    for(i in 1:p){
      for(j in 1:p){
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==1 | G.mat[j,i]==1)) TP <- TP + 1
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==0 | G.mat[j,i]==0)) FN <- FN + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==1 | G.mat[j,i]==1)) FP <- FP + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==0 | G.mat[j,i]==0)) TN <- TN + 1
      }
    }
  }
  prec <- TP / (TP + FP)
  if(TP+FP==0) prec <- 1
  
  rec <- TP / (TP + FN)
  if(TP+FN==0) rec <- 0
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  return(list(prec=prec, rec=rec, TPR=TPR, FPR=FPR))
}

