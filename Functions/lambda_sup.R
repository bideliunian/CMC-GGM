# A function that returns the minimum choice of lambda that guarantees full sparsity
lambda.sup <- function(data, p){
  n <- nrow(data)
  M <- ncol(data) / p
  candidates <- matrix(-1, p, p)
  for (i in 1:p) {
    Y <- data[,(i-1)*M + (1:M)]
    for(j in 1:p){
      if(j != i){
        X.Gj <- data[,(j-1)*M + (1:M)]
        candidates[i, j] <- norm(t(X.Gj)%*%Y,"F")/n 
      }
    }
  }
  return(max(candidates))
}
