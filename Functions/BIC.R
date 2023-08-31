########## BIC to choose tuning paramter ####################
bic <- function(S, Omega, M, n){
  blockwise_Frob <- function(x, M){
    # x: matrix
    # M: block size
    p <- dim(x)[1]
    if (p %% M != 0){
      stop("dimension of matrix cannot be divided by block size")
    }
    n.b <- p / M
    result <- matrix(0, nrow=n.b, ncol=n.b)
    for (i in 1:n.b){
      for (j in 1:n.b){
        result[i, j] <- norm(as.matrix(x[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]), "F")
      }
    }
    return(result)
  }
  
  Omega_block <- blockwise_Frob(Omega, M)
  diag(Omega_block) <- 0
  return(sum(diag(S%*%Omega)) - log(det(Omega)) + sum(Omega_block > 0) / 2 * M^2 *log(n))
}