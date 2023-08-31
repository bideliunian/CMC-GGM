# Extract XGj and get maximum 2-norm over j
max.2.norm <- function(X, M){
  n <- nrow(X)
  p <- ncol(X) / M
  norm.vec <- rep(0,p)
  for(j in 1:p){
    X.Gj <- X[,(j-1)*M + (1:M)]
    norm.vec[j] <- norm(X.Gj/sqrt(n), type="2")
  }
  return(max(norm.vec))
}

# Vectorize Y: from n*M to nm*1
Y.vectorize <- function(Y){
  n <- nrow(Y)
  M <- ncol(Y)
  res <- numeric(0)
  for(i in 1:n){
    res <- rbind(res, matrix(Y[i,],M,1))
  }
  res <- as.vector(res)
  return(res)
}

# Kronecker product X with I to get Z: nM*pM^2.
# Then combine a column of 1 as intercept
kprod <- function(X, M){
  Z <- kronecker(X, diag(M))
  return(Z)
}

# A.array to beta
col.vectorize <- function(A.array){
  R <- nrow(A.array)
  C <- ncol(A.array)
  res <- numeric(0)
  for(i in 1:C){
    res <- cbind(res, matrix(A.array[,i],R,1))
  }
  res <- as.vector(res)
  return(res)
}

# Iteration
group.lasso <- function(X, Y, lambda){
  # X and Y are both in matrix form. X: n*pM, Y: n*M
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X) / M
  
  Y.vec <- Y.vectorize(Y)
  Z <- kprod(X, M)
  Z.expand <- cbind(1,Z)
  grp.index <- c(NA, rep(1:p, each=M^2))
  
  fit <- grplasso(x=Z.expand, y=Y.vec, index=grp.index, lambda=lambda, model = LinReg())
  beta <- fit$coefficients
  epsilon <- Y.vec - Z.expand %*% beta
  epsilon.mat <- matrix(epsilon, nrow=n, ncol=M, byrow=T)
  beta.mat <- matrix(beta, nrow=, ncol=M, byrow=T)
  
  return(list(coef.mat = beta.mat, ))
}

